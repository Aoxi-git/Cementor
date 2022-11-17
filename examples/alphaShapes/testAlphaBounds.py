import random
N=20
for i in range(N):
    for j in range(N):
        for k in range(N):
            O.bodies.append(sphere((i+0.01*random.random(),j+0.01*random.random(),k+0.01*random.random()),0.5,color=(0.7,0.7,0.1)))
if False:
	O.bodies.append(sphere((9.50498,3485.33,9.50502),3465.82,color=(0.7,0.7,0.1)))

n = Vector3(1,1,1); n.normalize()
for b in O.bodies:
    if (b.state.pos).norm() < N/2:
        O.bodies.erase(b.id)
    if (b.state.pos - b.state.pos.dot(n)*n).norm()<2.5:
        O.bodies.erase(b.id)
        
TW=TesselationWrapper(far=100)
TW.triangulate()
a=TW.addBoundingPlane(1,True)
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

#from yade import qt
#qt.Controller()
#rr = qt.Renderer()
#graph = GlExtra_AlphaGraph()
#rr.extraDrawers = [graph]
#graph.alpha=4
#graph.shrinkedAlpha=((sqrt(graph.alpha)-0.5)**2)
#graph.fixedAlpha=True


v = qt.View()
