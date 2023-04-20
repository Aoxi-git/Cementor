from __future__ import print_function
from yade.beamconnectedsphere import *

O.engines = [
    ForceResetter(),
    InsertionSortCollider([Bo1_Sphere_Aabb()]),
    InteractionLoop(
        [Ig2_Sphere_Sphere_ScGeom()], 
        [Ip2_FrictMat_FrictMat_MindlinPhys(betan=0.0, betas=0.0)],
        [Law2_ScGeom_MindlinPhys_BeamConnectedSphere()]
    ),
    NewtonIntegrator(gravity=(0, 0, 0), label='newton', damping=0)
]

O.bodies.append(beamConnectedSphere([10, 10, 10], 2, wire=False, fixed=False))
O.bodies.append(beamConnectedSphere([13, 11, 10], 2, wire=False, fixed=False))
O.bodies[0].state.angVel = 5*Vector3(3, 1, 0)
O.bodies[1].state.angVel = 5*Vector3(3, 1, 0)
O.bodies[0].state.ori = Quaternion(0.894463, 0.291567, 0.172955, 0.291567)
O.bodies.append(beam(0, 1, damping=0.0))
O.dt = 1e-06
O.run(1000000, 1)  
A = O.bodies[0].state.pos - Vector3(10, 10, 10)
B = O.bodies[1].state.pos - Vector3(13, 11, 10)

if A.norm() != 0 or B.norm() != 0:
    raise YadeCheckError("Spinning Beam checktest: Law2_ScGeom_MindlinPhys_BeamConnectedSphere orientation calculation not correct.")