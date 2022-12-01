import random
N=10
for i in range(N):
    for j in range(N):
        for k in range(N):
            O.bodies.append(sphere((i+0.01*random.random(),j+0.01*random.random(),k+0.01*random.random()),0.49,color=(0.7,0.7,0.1)))

n = Vector3(1,1,1); n.normalize()
for b in O.bodies:
    if (b.state.pos).norm() < N/2:
        O.bodies.erase(b.id)
    if (b.state.pos - b.state.pos.dot(n)*n).norm()<2.5:
        O.bodies.erase(b.id)
        
boxes= O.bodies.append([box((N/2,-0.51,N/2),(N,0,N),fixed=True),box((N/2,N-0.49,N/2),(N,0,N),fixed=True)])

TW=TesselationWrapper(far=100)
TW.triangulate()
TW.addBoundingPlane(1,True)
TW.addBoundingPlane(1,False)
a = 2
shrinkedA = (sqrt(a)-0.5)**2


ag = TW.getAlphaGraph(alpha=a,shrinkedAlpha=shrinkedA,fixedAlpha=False)
graph = GlExtra_AlphaGraph(tesselationWrapper=TW,wire=False)
#graph.alpha=a
#graph.shrinkedAlpha=shrinkedA
#graph.fixedAlpha=False
graph.lighting=True

from yade import qt
rr = qt.Renderer()
rr.extraDrawers = [graph]

TW.applyAlphaForces(stress = -Matrix3.Identity, alpha=a,shrinkedAlpha=shrinkedA,fixedAlpha=False, reset = False)
for b in O.bodies: 
    if bool(b) and isinstance(b.shape,Sphere): b.shape.color *= O.forces.f(b.id).norm()
    else :  b.shape.color = (1,0,0)

#from yade import qt
#qt.Controller()
#rr = qt.Renderer()
#graph = GlExtra_AlphaGraph()
#rr.extraDrawers = [graph]
#graph.alpha=4
#graph.shrinkedAlpha=((sqrt(graph.alpha)-0.5)**2)
#graph.fixedAlpha=True

def updateMembrane():
    TW.triangulate(reset=True)
    TW.addBoundingPlane(1,True)
    TW.addBoundingPlane(1,False)
    TW.applyAlphaForces(stress = -1000*Matrix3.Identity, alpha=a,shrinkedAlpha=shrinkedA,fixedAlpha=False, reset = False)
    #FIXME: suboptimal, we can't re-use the above cosntruct because it would try to insert fictious spheres again
    # as a consequence we rebuild the whole thing just for display :-/
    TW.triangulate(reset=True)
    TW.addBoundingPlane(1,True)
    TW.addBoundingPlane(1,False)
    ag = TW.getAlphaGraph(alpha=a,shrinkedAlpha=shrinkedA,fixedAlpha=False)
    graph.refresh()


v = qt.View()


newton=NewtonIntegrator(damping=0.3)

O.engines=[
	ForceResetter(),
	InsertionSortCollider([Bo1_Sphere_Aabb(),Bo1_Box_Aabb()]),
	InteractionLoop(
		[Ig2_Sphere_Sphere_ScGeom(),Ig2_Box_Sphere_ScGeom()],
		[Ip2_FrictMat_FrictMat_FrictPhys()],
		[Law2_ScGeom_FrictPhys_CundallStrack()]
	),
	GlobalStiffnessTimeStepper(active=1,timeStepUpdateInterval=100,timestepSafetyCoefficient=0.8),
	#triax,
	newton,
	PyRunner(iterPeriod=50,command="updateMembrane()",label="membrane")
]
