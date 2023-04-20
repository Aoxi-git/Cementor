from yade.beamconnectedsphere import *
from yade import qt

# This example shows how to create a simply supported beam.
# Also illustrates how the beam-connected spheres can interact with other simulation objects.
# It shows how to add damping to the beam.
# It shows how to allow the beam to break according to Mohr Culomb's failure criteria.

# Set the simulation engine. Is important to define the engine Law before creating the beams. By default, when the
# beam is created it will create the interaction between the two spheres that are connected.
# If the engine is not defined, a warning saying that there's no law capable of handling the
# interaction will be shown. You can deactivate this behavior by setting the create_interaction=False
# when creating the beam. However, if your spheres are far apart, the interaction may not be created automatically. 
O.engines = [
        ForceResetter(),
        InsertionSortCollider([Bo1_Sphere_Aabb(), Bo1_Wall_Aabb()]),
        InteractionLoop(
                [Ig2_Sphere_Sphere_ScGeom(), Ig2_Wall_Sphere_ScGeom(label='ig2W')], 
                # This law works together with MindlinPhys, that's why we use the Ip2_FrictMat_FrictMat_MindlinPhys.
                [Ip2_FrictMat_FrictMat_MindlinPhys(en=0.75)],
                # When there's no beam connecting the spheres, the interaction is handled by Law2_ScGeom_MindlinPhys_Mindlin.
                [Law2_ScGeom_MindlinPhys_BeamConnectedSphere()] 
        ),
        NewtonIntegrator(gravity=(0, -9, 0), damping=0)
]

# Create materials.
O.materials.append(FrictMat(young=1e9, poisson=0.1, density=2500.0, frictionAngle=atan(0.4), label='Wmat'))
O.materials.append(FrictMat(young=1e6, poisson=0.1, density=2500.0, frictionAngle=atan(0.2), label='mat'))

# Create the beam-connected spheres.
L = 0.1  # Length of the beam.
n = 50   # Number of spheres -1 to use in the beam.
r = L/(2*n) # Radius of the spheres. We dont want the spheres to overlap because when the beam breaks, the overlap will create a big Hertz force.
color = [255. / 255., 102. / 255., 0. / 255.]

# Create beam-connected spheres.
nodeIds = [] # List of node ids so we can iterate over them later.
nodeIds.append(O.bodies.append(beamConnectedSphere([0, 0.1, 0], r, wire=False, fixed=True, material='mat', color=color)))
for i in range(1, n):
    nodeIds.append(O.bodies.append(beamConnectedSphere([i*L/n, 0.1, 0], r, wire=False, fixed=False, material='mat', color=color)))
nodeIds.append(O.bodies.append(beamConnectedSphere([L, 0.1, 0], r, wire=False, fixed=True, material='mat', color=color)))

# Connect the spheres with beams. 
# The numbers are the indices of the spheres in the O.bodies list to be connected.
# Only beam-connected spheres can be connected with beams.
for i, j in zip(nodeIds[:-1], nodeIds[1:]): # Iterate over the list of node ids by pairs
    O.bodies.append(beam(i, j, damping=0.05, fracture=True, phi=0.2, cohesion=6.2e4)) #Create the beams. We use Rayleigh damping and Mohr-Coulomb's failure criteria.

# Add a sphere that will collide with the beam and break it.
# Beam-connected spheres can interact with spheres.
O.bodies.append(utils.sphere([L/2, 0.1 + L, 0.0], 3.2*r, material='mat'))

# Add a wall
# Beam-connected spheres can interact with walls.
O.bodies.append(utils.wall([0,0,0], 1, material=O.materials["Wmat"]))

# Render stripes in the spheres so we can see them rotate.
Gl1_Sphere.stripes=1

# Set the time step.
O.dt = 0.1*PWaveTimeStep()
O.stopAtIter = 270000

# Render the scene.
O.saveTmp()
yade.qt.Renderer()
qt.View()