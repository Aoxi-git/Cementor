from yade.beamconnectedsphere import *
from yade import qt
import yade.math as m
import yade.minieigenHP as me

# This example shows how to create beam between neigoring spheres.
# Also illustrates how the beam-connected spheres
# can interact with other simulation objects.
# It shows how to add damping to the beam.
# It shows how to allow the beam to break according to Mohr Culomb's failure criteria.
# It shows how to use YADE high precision functions with the beams connected spheres.

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
                [Ip2_FrictMat_FrictMat_MindlinPhys(en=m.toHP1('0.75'))],
                # When there's no beam connecting the spheres, the interaction is handled by Law2_ScGeom_MindlinPhys_Mindlin.
                [Law2_ScGeom_MindlinPhys_BeamConnectedSphere()] 
        ),
        NewtonIntegrator(gravity=me.HP1.Vector3(0, -9, 0), damping=m.toHP1('0'))
]

# Create materials.
O.materials.append(FrictMat(young=m.toHP1('1e40'), poisson=m.toHP1('0.1'), density=m.toHP1('2500.0'), frictionAngle=m.atan(m.toHP1('0.4')), label='Wmat'))
O.materials.append(FrictMat(young=m.toHP1('5e6'), poisson=m.toHP1('0.1'), density=m.toHP1('500.0'), frictionAngle=m.atan(m.toHP1('0.2')), label='mat'))


# Create a quaternion that rotates the spheres 0.5 rad around all the axis.
#                            w,        i,        j,       k
q = me.HP1.Quaternion(0.894463, 0.291567, 0.172955, 0.291567)

# Create beam connected spheres arranged in a cube
nx = 6
ny = 6
nz = 6
r = m.toHP1('0.1')
for i in range(nx):
    for j in range(ny):
        for k in range(nz):
            O.bodies.append(beamConnectedSphere(q*me.HP1.Vector3(i*2*r, j*2*r, k*2*r), radius=r, material='mat'))

# Create beams between neighboring spheres.
buffer = r/5

# Checks if two spheres are in contact. If they are, it creates a beam between them.
# Also, if adding the buffer to the radius of the spheres, they are still in contact, it creates a beam between them.
# We allow the beam to break according to Mohr Culomb's failure criteria.
# We add Rayleigh damping to the beam.
BeamConnecNeighborSpheres(buffer, damping=0.05, fracture=True, phi=m.toHP1('0.1'), cohesion=m.toHP1('1e5'))

# Add a wall
# Beam-connected spheres can interact with walls.
O.bodies.append(utils.wall(me.HP1.Vector3(0,-1.0,0), 1, material=O.materials["Wmat"]))

# Render stripes in the spheres so we can see them rotate.
Gl1_Sphere.stripes=1

# Set the time step.
O.dt = 8e-5
#O.stopAtIter = 270000

# Render the scene.
O.saveTmp()
yade.qt.Renderer()
v = qt.View()
v.sceneRadius = 5