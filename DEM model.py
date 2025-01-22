import pandas as pd
import numpy as np
from yade import pack, geom, qt, plot, utils
from yade import FrictMat, ForceResetter, InsertionSortCollider, InteractionLoop, NewtonIntegrator, PyRunner
from yade.utils import radians

# Load the dataset
file_path = r'C:\Users\ENG. JMK\OneDrive\Desktop\GEOCHEMICAL DATA\DataSet_Thaba_Classification.csv'
geochemical_data = pd.read_csv(file_path, delimiter=";")

# Filter relevant columns
relevant_columns = ["DepthFrom", "DepthTo", "SiO2_%", "FeO_%", "MgO_%", "Al2O3_%", "CaO_%", "Stratigraphy"]
filtered_data = geochemical_data[relevant_columns]

# Calculate material properties (example: density)
filtered_data["SiO2_%"] = pd.to_numeric(filtered_data["SiO2_%"], errors="coerce")
filtered_data["Density"] = 2.7 * (1 + 0.01 * filtered_data["SiO2_%"])

# YADE: Define slope geometry and materials
slope_material = FrictMat(young=1e9, poisson=0.3, density=2700, frictionAngle=radians(30))
O.materials.append(slope_material)

# Create a box for the slope
O.bodies.append(geom.facetBox(center=(50, 50, 50), extents=(50, 50, 50), wallMask=31, material=slope_material))

# Create particles
sp = pack.SpherePack()
sp.makeCloud(minCorner=(0, 0, 0), maxCorner=(100, 100, 100), rMean=1, rRelFuzz=0.3)
sp.toSimulation(material=slope_material)

# Define engines
def addData():
    if O.bodies:  # Ensure there are bodies in the simulation
        displacement = O.bodies[0].state.pos[2] - O.bodies[0].state.refPos[2]
        plot.addData(t=O.time, displacement=displacement)

plot.plots = {'t': ('displacement',)}
O.engines = [
    ForceResetter(),
    InsertionSortCollider([Bo1_Sphere_Aabb(), Bo1_Facet_Aabb()]),
    InteractionLoop(
        [Ig2_Sphere_Sphere_ScGeom(), Ig2_Facet_Sphere_ScGeom()],
        [Ip2_FrictMat_FrictMat_FrictPhys()],
        [Law2_ScGeom_FrictPhys_CundallStrack()]
    ),
    NewtonIntegrator(damping=0.2, gravity=(0, 0, -9.81)),
    PyRunner(command="addData()", iterPeriod=100)
]

# Run the simulation
O.run(10000, wait=True)
plot.plot()
plt.show()
