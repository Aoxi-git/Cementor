import yade.math as mth
import yade.minieigenHP as mne

Body = utils.sphere(mne.Vector3(0, 0, 0), 1.0)
ID = O.bodies.append(Body)
O.bodies[ID].state.mass = mth.Real("1.0")
O.bodies[ID].state.inertia = mne.Vector3(1.0, 1.0, 1.0)

O.engines=[
    ForceResetter(),
    InsertionSortCollider(),
    InteractionLoop(),
    NewtonIntegrator(damping=0.0, gravity=mne.Vector3(0, 0, 0)),
]

O.forces.setPermF(ID, mne.Vector3(1.0, 0, 0))
O.forces.setPermT(ID, mne.Vector3(1.0, 0, 0))

O.stopAtTime = mth.Real("1.0")
O.dt = mth.Real("1.0e-4")

O.run()
O.wait()

v = O.bodies[ID].state.vel[0]
w = O.bodies[ID].state.angVel[0]

if abs(v - mth.Real("1.0001")) > mth.Real("1e-10"):
    raise YadeCheckError("setPermF not working, expected vel = 1.0001, got vel = ", v)

if abs(w - mth.Real("1.0001")) > mth.Real("1e-10"):
    raise YadeCheckError("setPermT not working, expected angVel = 1.0001, got angVel = ", w)
