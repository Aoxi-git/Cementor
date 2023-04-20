from yade.beamconnectedsphere import *
from yade import qt

# This example shows how to create a beam in a cantilever.
# Also illustrates how the beam-connected spheres
# can interact with other simulation objects.
# It shows how to add damping to the beam.
# It shows that the beam connection supports periodic boundaries.

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
        NewtonIntegrator(gravity=(0, -9, 0), damping=0),
]

# Create a material.
O.materials.append(FrictMat(young=4e8, poisson=0.1, density=2500.0, frictionAngle=atan(0.2), label='mat'))

# Create periodic boundaries.
O.periodic = True
O.cell.hSize = Matrix3(0.14, 0, 0, 0, 0.18, 0, 0, 0, 0.06)
O.cell.velGrad = Matrix3(-.1, 0, 0, 0, -.1, 0, 0, 0, -.1)

# Create the beam-connected spheres.
L = 0.1  # Length of the beam.
n = 50   # Number of spheres to use in the beam.
r = L /(1.1*n) # Radius of the spheres.
color = [255. / 255., 102. / 255., 0. / 255.]

# Create beam-connected spheres.
O.bodies.append(beamConnectedSphere([0.07, 0.05, 0.05], r, wire=False, fixed=True, material='mat', color=color))
for i in range(1, n):
    O.bodies.append(beamConnectedSphere([0.07 + i*L/n, 0.05, 0.05], r, wire=False, fixed=False, material='mat', color=color))

# Connect the spheres with beams. 
BeamConnecNeighborSpheres(0, damping=0.05)

# To spice things up let's add a sphere that will collide with the beam.
# Beam-connected spheres can interact with spheres.
O.bodies.append(utils.sphere([0.07 + (n-3)*L/n, 0.05 + L, 0.05], 3*r, material='mat'))

# Render stripes in the spheres so we can see them rotate.
Gl1_Sphere.stripes=1

# Set the time step.
O.dt = 0.5*PWaveTimeStep()
O.stopAtIter = 368000

# Render the scene.
O.saveTmp()
yade.qt.Renderer()
qt.View()