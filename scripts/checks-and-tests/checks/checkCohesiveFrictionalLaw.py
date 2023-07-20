from numpy import array
#from yade import plot   ## uncomment plotting lines for debugging

O.periodic = True
O.cell.hSize = Matrix3(10, 0, 0, 0, 10, 0, 0, 0, 10)

results=[]
def getState():
	state = [int(O.interactions.has(0,1) and O.interactions[0,1].isReal)]
	if state[0]:
		i = O.interactions[0,1]
		for x in [int(i.phys.cohesionBroken),i.phys.normalForce.norm(),i.phys.shearForce.norm(),i.phys.moment_bending.norm(),i.phys.moment_twist.norm()]:
			state.append(x)
	results.append(numpy.array(state))
	#print(state)

expectedResults = [array([1., 0., 0., 0., 0., 0.]),
 array([  1. ,   0. , 399.6,   0. ,   0. ,   0. ]),
 array([  1.,   0., 999.,   0.,   0.,   0.]),
 array([False]),
 array([False]),
 array([ 1. ,  1. , 98.9,  0. ,  0. ,  0. ]),
 array([  1. ,   1. , 198.8,   0. ,  99.9,   0. ]),
 array([  1. ,   1. , 198.8,   0. , 198.8,   0. ]),
 array([  1. ,   1. , 198.8,   0. , 198.8,   0. ]),
 array([  1. ,   1. , 198.8,   0. ,  1. ,   0. ]),
 array([  1. ,   1. , 198.8,   0. ,  1.00497762, 99.9]),
 array([  1. ,   1. , 198.8,   0. ,  1.00998001, 198.8]),
 array([  1. ,   1. , 198.8,   0. ,  1.0150073, 198.8]),
 array([  1. ,   0. ,   0. ,   0. , 249.75,   0. ]),
 array([False]),
 array([  1.,    0.,  999.  ,   0.,   0.,   0.]),
 array([1.0000e+00, 0.0000e+00, 1.0989e+03, 0.0000e+00, 1.2987e+03, 0.0000e+00]),
 array([1.0000e+00, 1.0000e+00, 1.0989e+03, 0.0000e+00, 1.0989e+03, 0.0000e+00])]


O.engines = [  # define engines, main functions for simulation
	ForceResetter(),
	InsertionSortCollider([Bo1_Sphere_Aabb()],
						label='ISCollider'),
	InteractionLoop(
		[Ig2_Sphere_Sphere_ScGeom6D()],
		[Ip2_CohFrictMat_CohFrictMat_CohFrictPhys(setCohesionOnNewContacts=True,label="ip2")],
		[Law2_ScGeom6D_CohFrictPhys_CohesionMoment(useIncrementalForm=True,always_use_moment_law=True)]
	),
	#PyRunner(command='plot.addData(ii=O.iter,x=O.bodies[1].state.pos[0]-7,Fx=O.forces.f(0)[0],Mz=O.forces.t(0)[2])', iterPeriod=1,initRun=True),
	NewtonIntegrator(label='newton')
]

O.dt=1
O.materials.append( CohFrictMat(young=1e3, normalCohesion=1e3, shearCohesion=0.5e3, etaRoll=1, etaTwist=1, momentRotationLaw=True, fragile=True) )
O.bodies.append( sphere((4,5,5),1,fixed=True))
O.bodies.append( sphere((7,5,5),1,fixed=True))

i = createInteraction(0, 1)
ip2.setCohesionOnNewContacts = False
i.phys.unp = -(O.bodies[i.id1].state.pos - O.bodies[i.id2].state.pos).norm() + O.bodies[i.id2].shape.radius + O.bodies[i.id1].shape.radius

O.saveTmp()
#plot.plots={'ii':('Fx',None,'Mz')}
#plot.plot()

O.reload()
print("Pull contact until it breaks")
getState()
O.bodies[1].state.vel=(0.0999,0,0)
O.run(5,True)
getState()
O.run(6,True)
getState()
O.run(1,True)
getState()
print("Inverse velocity to create new non-cohesive contact")
O.bodies[1].state.vel= -O.bodies[1].state.vel
O.run(23,True)
getState()
O.run(1,True)
getState()
print("Apply pure bending")
O.bodies[1].state.vel=(0,0,0)
O.bodies[1].state.angVel=(0,0,0.0999)
O.bodies[0].state.angVel=(0,0,-0.0999)
O.run(1,True)
getState()
O.run(5,True)
getState()
O.run(1,True)
getState()
print("Un-bend")
O.bodies[1].state.angVel=(0,0,-0.0999)
O.bodies[0].state.angVel=(0,0,0.0999)
O.run(2,True)
getState()
print("Twist")
O.bodies[1].state.angVel=(0,0,0)
O.bodies[0].state.angVel=(2*0.0999,0,0)
O.run(1,True)
getState()
O.run(1,True)
getState()
O.run(1,True)
getState()
print("Break remote interaction by bending")
O.reload()
O.bodies[1].state.angVel=(0,0,-0.00999)
O.bodies[0].state.angVel=(0,0,0.00999)
O.run(25,True)
getState()
O.run(1,True)
getState()
print("Bend compressed contact with residual friction")
O.reload()
O.bodies[1].state.vel=(-0.0999,0,0)
O.run(11,True)
getState()
O.bodies[1].state.vel=(0,0,0)
O.bodies[1].state.angVel=(0,0,-0.0999)
O.bodies[0].state.angVel=(0,0,0.0999)
O.run(13,True)
getState()
O.run(1,True)
getState()

difference = numpy.linalg.norm((numpy.linalg.norm(numpy.array(results) - numpy.array(expectedResults))))
if difference>1e-8: 
	raise YadeCheckError('ThermalEngine checktest: fluid temp incorrect (difference='+str(difference)+')')
else:
	print("CohesiveFrictional model passed with difference="+str(difference))
