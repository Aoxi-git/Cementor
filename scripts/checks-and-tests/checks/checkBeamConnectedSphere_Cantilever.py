from yade.beamconnectedsphere import *
import yade.math as m

def solution(x):
    return a*x**4 + b*x**3 + c*x**2

# SI units
pi = m.Pi()
g = -9.81
n = 80
L = 0.1
r = L/(2.0*n)
E = 1.0e8
D = 800.0
m = D*(4/3)*pi*r**3
A = pi*r**2
I = A*A/(4.0*pi)     

w = m*g*n/L

a = w/(24*E*I)
b = - 24*a*L/6
c = -(6*b*L+12*a*L*L)/2

O.engines = [
    ForceResetter(),
    InsertionSortCollider([Bo1_Sphere_Aabb()]),
    InteractionLoop(
        [Ig2_Sphere_Sphere_ScGeom()], 
        [Ip2_FrictMat_FrictMat_MindlinPhys(betan=0.0, betas=0.0)],
        [Law2_ScGeom_MindlinPhys_BeamConnectedSphere()] 
    ),
    NewtonIntegrator(gravity=(0, g, 0), damping=0)
]

O.materials.append(FrictMat(young=E, density=D, label='mat'))

O.bodies.append(beamConnectedSphere([0, 0, 0], r, fixed=True, material='mat'))
for i in range(1, n):
    O.bodies.append(beamConnectedSphere([i*L/n, 0, 0], r, material='mat'))

BeamConnecNeighborSpheres(r/3, damping=0.35)

O.dt = 1e-06
O.run(400000, 1)

error = []
for i in O.bodies:
    if isinstance(i.shape, BeamConnectedSphere):
        x = i.state.pos[0]
        y = i.state.pos[1]
        err = (y - solution(x))**2
        error.append(err)

Error = sum(error)/len(error)

if Error > 1e-8:
    raise YadeCheckError("Cantilever Beam checktest: Law2_ScGeom_MindlinPhys_BeamConnectedSphere Force model incorrect.")