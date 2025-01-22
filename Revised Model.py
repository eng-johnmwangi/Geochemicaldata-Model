import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from yade import plot, utils, pack, geom, qt, FrictMat
from math import radians

# Step 1: Load and Preprocess Data
data = pd.read_csv("DataSet_Thaba_Classification.csv", delimiter=";")

# Filter relevant columns
relevant_columns = ["DepthFrom", "DepthTo", "SiO2_%", "FeO_%", "MgO_%", "Al2O3_%", "CaO_%", "Stratigraphy"]
filtered_data = data[relevant_columns]

# Calculate material properties (example: density)
filtered_data["Density"] = 2.7 * (1 + 0.01 * filtered_data["SiO2_%"])  # Example formula

# Step 2: Define Slope Geometry
# Create a simple 2D slope geometry
slope_height = 100  # meters
slope_angle = 30  # degrees

# Step 3: Set Up DEM Simulation in Yade
O.materials.append(FrictMat(young=1e9, poisson=0.3, density=2700, frictionAngle=radians(30)))

# Create a box for the slope
O.bodies.append(geom.facetBox((0, 0, 0), (100, 100, 100)))

# Create particles for the slope material
sp = pack.SpherePack()
sp.makeCloud(minCorner=(0, 0, 0), maxCorner=(100, 100, 100), rMean=1, rRelFuzz=0.3)
sp.toSimulation()

# Initialize reference positions
for b in O.bodies:
    b.state.refPos = b.state.pos

# Step 4: Incorporate Pore Pressure Effects
def calculate_pore_pressure(z, P0=0, gamma_w=9.81):
    """
    Calculate pore pressure at depth z.
    P0: Surface pore pressure (default = 0).
    gamma_w: Unit weight of water (default = 9.81 kN/m^3).
    """
    return P0 + gamma_w * z

def update_effective_stress(bodies, P0=0, gamma_w=9.81):
    """
    Update effective stress for all particles based on pore pressure.
    """
    for b in bodies:
        z = b.state.pos[2]  # Depth of the particle
        pore_pressure = calculate_pore_pressure(z, P0, gamma_w)
        b.state.customDict['effective_stress'] = b.state.stress - pore_pressure

# Step 5: Run Simulation
O.engines = [
    ForceResetter(),
    InsertionSortCollider([Bo1_Sphere_Aabb(), Bo1_Facet_Aabb()]),
    InteractionLoop(
        [Ig2_Sphere_Sphere_ScGeom(), Ig2_Facet_Sphere_ScGeom()],
        [Ip2_FrictMat_FrictMat_FrictPhys()],
        [Law2_ScGeom_FrictPhys_CundallStrack()]
    ),
    NewtonIntegrator(damping=0.2, gravity=(0, 0, -9.81)),
    PyRunner(command="addData()", iterPeriod=100),
    PyRunner(command="update_effective_stress(O.bodies)", iterPeriod=100)
]

# Function to collect data during simulation
def addData():
    plot.addData(
        t=O.time,
        displacement=O.bodies[0].state.pos[2] - O.bodies[0].state.refPos[2]
    )

# Step 6: Run Simulation and Plot Results
O.run(10000, wait=True)
plt.figure()
plt.plot(plot.data['t'], plot.data['displacement'])
plt.xlabel("Time (s)")
plt.ylabel("Displacement (m)")
plt.title("Displacement vs. Time")
plt.show()
