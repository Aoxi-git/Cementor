# -*- encoding=utf-8 -*-


# Abstract : this script defines a "rock joint" sample : two rectangular blocks separated by an horizontal joint surface. Imposing relative movements of the blocks, that are clumps, allows to test directly the behaviour of the joint, described by JCFpm model.
# jerome.duriez@3sr-grenoble.fr


# Mechanical properties of rock matrix and rock joint :
from __future__ import print_function
from past.builtins import execfile
def mat(): return JCFpmMat(type=1,young=15.e9,frictionAngle=radians(35),density=3000,poisson=0.35,tensileStrength=4.5e6,cohesion=45.e6,jointNormalStiffness=5.e7,jointShearStiffness=2.5e7,jointCohesion=0.,jointTensileStrength=0.,jointFrictionAngle=radians(35.),jointDilationAngle=0.0)


# --- Creating a sample of spheres

# definition of a predicate 
from yade import pack
Lx = 10
Ly = 10
Lz = 6
pred = pack.inAlignedBox((0,0,0),(Lx,Ly,Lz))
# use of randomDensePack() function
nSpheres = 1500.0
poros=0.13 # apparently the value of porosity of samples generated by pack.randomDensePack
rMeanSpheres = pow(Lx*Ly*Lz*3.0/4.0*(1-poros)/(pi*nSpheres),1.0/3.0)
print('\nGenerating sphere sample, be patient')
sp = pack.randomDensePack(pred,radius=rMeanSpheres,rRelFuzz=0.3,memoizeDb='/tmp/gts-triax-packings.sqlite',returnSpherePack=True)
sp.toSimulation(color=(0.9,0.8,0.6),wire=False,material=mat)
print('Sphere sample generated !')


# --- The joint surface : half of the height
import gts
v1 = gts.Vertex(0 , 0 , Lz/2.0)
v2 = gts.Vertex(Lx, 0 , Lz/2.0)
v3 = gts.Vertex(Lx, Ly, Lz/2.0)
v4 = gts.Vertex(0 , Ly, Lz/2.0)

e1 = gts.Edge(v1,v2)
e2 = gts.Edge(v2,v4)
e3 = gts.Edge(v4,v1)
f1 = gts.Face(e1,e2,e3)

e4 = gts.Edge(v4,v3)
e5 = gts.Edge(v3,v2)
f2 = gts.Face(e2,e4,e5)

s1 = gts.Surface()
s1.add(f1)
s1.add(f2)

facet = gtsSurface2Facets(s1,wire = False,material=mat)
O.bodies.append(facet)


# --- Identification of spheres onJoint, and so on:
execfile('identifBis.py')


# --- Engines definition
O.engines=[
	ForceResetter(),
	InsertionSortCollider([Bo1_Sphere_Aabb()]),
	InteractionLoop(
		[Ig2_Sphere_Sphere_ScGeom()],
		[Ip2_JCFpmMat_JCFpmMat_JCFpmPhys(cohesiveTresholdIteration=1)],
		[Law2_ScGeom_JCFpmPhys_JointedCohesiveFrictionalPM(smoothJoint=True)]),
	GlobalStiffnessTimeStepper(timestepSafetyCoefficient=0.8),
        NewtonIntegrator(damping=0.2),
        PyRunner(command='afficheIt()',initRun=True,iterPeriod=1000),
]
def afficheIt():
	print('It', O.iter)

O.step()

# --- Clumping the blocks
upperBlock=[]
lowerBlock=[]
for inte in O.interactions:
    if not inte.phys.isOnJoint:
       bod1 = O.bodies[inte.id1]
       bod2 = O.bodies[inte.id2]
       if bod1.state.pos[2]<Lz/2.0:
          if not (bod1.id in lowerBlock):
                  lowerBlock.append(bod1.id)
                  bod1.shape.color=Vector3(1,0,0)
          if bod2.state.pos[2]>Lz/2.0:
             print('\n **** ERROR !!!! ******* \n\n')
          else:
             if not (bod2.id in lowerBlock):
                     lowerBlock.append(bod2.id)
                     bod2.shape.color=Vector3(1,0,0)
       else:
          if not (bod1.id in upperBlock):
                  upperBlock.append(bod1.id)
                  bod1.shape.color=Vector3(0,0,1)
          if bod2.state.pos[2]<Lz/2.0:
             print('\n **** ERROR !!!! ******* \n\n')
          else:
                  if not (bod2.id in upperBlock):
                          upperBlock.append(bod2.id)
                          bod2.shape.color=Vector3(0,0,1)
print('\n Clumping upper block, be patient')
idUpperClump=O.bodies.clump(upperBlock)
print('Clumped !')

print('\n Clumping lower block, be patient')
idLowerClump=O.bodies.clump(lowerBlock)
print('Clumped !')

upperClump = O.bodies[idUpperClump]
lowerClump = O.bodies[idLowerClump]
z0 = upperClump.state.pos[2]
upperClump.dynamic=False
lowerClump.dynamic=False


# --- Saving data
O.engines =O.engines+[PyRunner(command='dataCollector()',initRun=True,iterPeriod=500)]
def dataCollector():
	plot.addData(iterations=O.iter,un=-(upperClump.state.pos[2]-z0),Fz=O.forces.f(idUpperClump,sync=True)[2],Ft=sqrt(pow(O.forces.f(idUpperClump,sync=True)[1],2.)+pow(O.forces.f(idUpperClump,sync=True)[0],2.)),Fzinf=O.forces.f(idLowerClump,sync=True)[2],uf=unbalancedForce(),ec=kineticEnergy())

from yade import plot
plot.plots={'un':('Fz')}

# --- Simulation !
compSpeed = 0.0075 # put here a positive value, for compression
upperClump.state.vel=Vector3(0,0,-compSpeed)
nIt = 10000
print('\nComputation begins now, for', nIt,'iterations')
O.run(nIt,wait=True)
print('Computation just finished ! \n')

yade.qt.View()
yade.qt.Controller()
plot.plot()
