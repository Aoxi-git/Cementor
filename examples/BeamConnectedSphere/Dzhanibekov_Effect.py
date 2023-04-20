from yade.beamconnectedsphere import *
from yade import qt

# This example shows the Dzhanibekov effect.
# Here we show how to create beam-connected spheres,
# how to connect them with beams, and how to run the simulation.

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
        NewtonIntegrator(gravity=(0, 0, 0), damping=0)
]

# Create a material.
O.materials.append(FrictMat(young=1e9, poisson=0.1, density=1000.0, label='mat'))

# Create the beam-connected spheres. They work like spheres. But, they have some extra methods to handle the connections.
O.bodies.append(beamConnectedSphere([10, 10, 10], 2, material='mat'))
O.bodies.append(beamConnectedSphere([13.8, 10, 10], 2, material='mat'))
O.bodies.append(beamConnectedSphere([17.6, 10, 10], 2, material='mat'))
O.bodies.append(beamConnectedSphere([13.8, 10, 13.8], 2, material='mat'))

# Set the initial velocities, make the spheres spin.
v = 66.6
O.bodies[0].state.vel = Vector3(0, -v, 0)
O.bodies[2].state.vel = Vector3(0, v, 0)
O.bodies[1].state.angVel = Vector3(0.1, 0, v/4)
O.bodies[3].state.angVel = Vector3(0.1, 0, v/4)

# Connect the spheres with beams. 
# The numbers are the indices of the spheres in the O.bodies list to be connected.
# Only beam-connected spheres can be connected with beams.
O.bodies.append(beam(0, 1))
O.bodies.append(beam(1, 2))
O.bodies.append(beam(1, 3))

# Render stripes in the spheres so we can see them rotate.
Gl1_Sphere.stripes=1

# Set the time step.
O.dt = 3e-6

# Render the scene.
O.saveTmp()
yade.qt.Renderer()
qt.View()