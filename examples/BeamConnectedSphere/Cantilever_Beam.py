from yade.beamconnectedsphere import *
from yade import qt

# This example shows how to create a beam in a cantilever.
# Also illustrates how the beam-connected spheres
# can interact with other simulation objects.
# It shows how to add damping to the beam.

# Set the simulation engine. Is important to define the engine Law before creating the beams. By default, when the
# beam is created it will create the interaction between the two spheres that are connected.
# If the engine is not defined, a warning saying that there's no law capable of handling the
# interaction will be shown. You can deactivate this behavior by setting the create_interaction=False
# when creating the beam. However, if your spheres are far apart, the interaction may not be created automatically. 
O.engines = [
        ForceResetter(),
        InsertionSortCollider([Bo1_Sphere_Aabb()]),
        InteractionLoop(
                [Ig2_Sphere_Sphere_ScGeom()], 
                # This law works together with MindlinPhys, that's why we use the Ip2_FrictMat_FrictMat_MindlinPhys.
                [Ip2_FrictMat_FrictMat_MindlinPhys(betan=0.0, betas=0.0)],
                # When there's no beam connecting the spheres, the interaction is handled by Law2_ScGeom_MindlinPhys_Mindlin.
                [Law2_ScGeom_MindlinPhys_BeamConnectedSphere()] 
        ),
        NewtonIntegrator(gravity=(0, -9, 0), damping=0)
]

# Create a material.
O.materials.append(FrictMat(young=4e8, poisson=0.1, density=2500.0, frictionAngle=atan(0.2), label='mat'))

# Create the beam-connected spheres.
L = 0.1  # Length of the beam.
n = 50   # Number of spheres to use in the beam.
r = L /n # Radius of the spheres.
color = [255. / 255., 102. / 255., 0. / 255.]

# Create beam-connected spheres.
nodeIds = [] # List of node ids so we can iterate over them later.
nodeIds.append(O.bodies.append(beamConnectedSphere([0, 0, 0], r, wire=False, fixed=True, material='mat', color=color)))
for i in range(1, n):
    nodeIds.append(O.bodies.append(beamConnectedSphere([i*L/n, 0, 0], r, wire=False, fixed=False, material='mat', color=color)))

# Connect the spheres with beams. 
# The numbers are the indices of the spheres in the O.bodies list to be connected.
# Only beam-connected spheres can be connected with beams.
for i, j in zip(nodeIds[:-1], nodeIds[1:]): # Iterate over the list of node ids by pairs
    O.bodies.append(beam(i, j, damping=0.05)) #Create the beams. We use Rayleigh damping

# To spice things up let's add a sphere that will collide with the beam.
# Beam-connected spheres can interact with spheres.
O.bodies.append(utils.sphere([L/2, L, 0.0], 3*r, material='mat'))

# Render stripes in the spheres so we can see them rotate.
Gl1_Sphere.stripes=1

# Set the time step.
O.dt = 0.4*PWaveTimeStep()
O.stopAtIter = 751000

# Render the scene.
O.saveTmp()
yade.qt.Renderer()
qt.View()