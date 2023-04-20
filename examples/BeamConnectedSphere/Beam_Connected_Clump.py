from yade.beamconnectedsphere import *
from yade import qt

# This example shows how to create clumbs of beam-connected spheres.
# Also illustrates how to connect the clumps with beams.
# It shows how to add damping to the beam.
# It shows how to allow the beam to break according to Mohr Culomb's failure criteria.
# It shows that the beam-connected spheres can interact with facets.

# Set the simulation engine. Is important to define the engine Law before creating the beams. By default, when the
# beam is created it will create the interaction between the two spheres that are connected.
# If the engine is not defined, a warning saying that there's no law capable of handling the
# interaction will be shown. You can deactivate this behavior by setting the create_interaction=False
# when creating the beam. However, if your spheres are far apart, the interaction may not be created automatically. 
O.engines = [
        ForceResetter(),
        InsertionSortCollider([Bo1_Sphere_Aabb(), Bo1_Facet_Aabb()]),
        InteractionLoop(
                [Ig2_Sphere_Sphere_ScGeom(), Ig2_Facet_Sphere_ScGeom(label='ig2W')], 
                # This law works together with MindlinPhys, that's why we use the Ip2_FrictMat_FrictMat_MindlinPhys.
                [Ip2_FrictMat_FrictMat_MindlinPhys(en=0.75)],
                # When there's no beam connecting the spheres, the interaction is handled by Law2_ScGeom_MindlinPhys_Mindlin.
                [Law2_ScGeom_MindlinPhys_BeamConnectedSphere()] 
        ),
        NewtonIntegrator(gravity=(0, -9, 0), damping=0)
]

# Create materials.
O.materials.append(FrictMat(young=1e7, poisson=0.1, density=500.0, frictionAngle=atan(0.2), label='mat'))
O.materials.append(FrictMat(young=1e9, poisson=0.1, density=2500.0, frictionAngle=atan(0.2), label='Wmat'))

# Create the 2 clumps of beam-connected spheres.
r = 0.1 # Radius of the spheres.
H = 2   # Height of the clumps.
dz = 10.5*r # Distance between the clumps.

O.bodies.appendClumped(
    [beamConnectedSphere([0, H, 0], material=O.materials['mat'], radius=3*r), 
    beamConnectedSphere([2*r, H, 0], material=O.materials['mat'], radius=3*r),
    beamConnectedSphere([5*r, H, 0], material=O.materials['mat'], radius=3*r)])

O.bodies.appendClumped(
    [beamConnectedSphere([dz, H, 0], material=O.materials['mat'], radius=3*r), 
    beamConnectedSphere([dz + 2*r, H + 3*r, 0], material=O.materials['mat'], radius=3*r),
    beamConnectedSphere([dz + 5*r, H + 5*r, 0], material=O.materials['mat'], radius=3*r)])

# Create the beam connecting the clumps.
# The numbers are the indices of the spheres in the O.bodies list to be connected.
# Only beam-connected spheres can be connected with beams.
O.bodies.append(beam(2, 4, damping=0.05, fracture=False, phi=0.3, cohesion=1e5)) #Create the beams. We use Rayleigh damping and Mohr-Coulomb's failure criteria.

# Create the facets to work as the flor.
L = 9
O.bodies.append(
    [
    facet([(-L, 0.0, -L), (L, 0.0, L), (L, 0.0, -L)], material='Wmat'),
    facet([(L, 0.0, L), (-L, 0.0, -L), (-L, 0.0, L)], material='Wmat'),
    ])

# Render stripes in the spheres so we can see them rotate.
Gl1_Sphere.stripes=1

# Set the time step.
O.dt = 1e-6
O.stopAtIter = 4510000

# Render the scene.
O.saveTmp()
yade.qt.Renderer()
qt.View()