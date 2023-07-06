if (('PFVFLOW' in features) and ('TWOPHASEFLOW' in features)):
	from yade import pack
	from yade import export
	from yade import timing
	from yade import plot
	import time
	from math import *
	import builtins

	num_spheres = 1000  # number of spheres
	young = 1e6
	compFricDegree = 3  # initial contact friction during the confining phase
	finalFricDegree = 30  # contact friction during the deviatoric loading
	mn, mx = Vector3(0, 0, 0), Vector3(1, 1, 0.4)  # corners of the initial packing
	graindensity = 2600
	toleranceWarning = 1.e-11
	toleranceCritical = 1.e-6

	O.materials.append(FrictMat(young=young, poisson=0.5, frictionAngle=radians(compFricDegree), density=graindensity, label='spheres'))
	O.materials.append(FrictMat(young=young, poisson=0.5, frictionAngle=0, density=0, label='walls'))
	walls = aabbWalls([mn, mx], thickness=0, material='walls')
	wallIds = O.bodies.append(walls)

	sp = pack.SpherePack()
	sp.makeCloud(mn, mx, -1, 0.3333, num_spheres, False, 0.95, seed=1)  #"seed" make the "random" generation always the same
	sp.toSimulation(material='spheres')

	triax = TriaxialStressController(
	        maxMultiplier=1. + 2e4 / young,  # spheres growing factor (fast growth)
	        finalMaxMultiplier=1. + 2e3 / young,  # spheres growing factor (slow growth)
	        thickness=0,
	        stressMask=7,
	        max_vel=0.005,
	        internalCompaction=True,  # If true the confining pressure is generated by growing particles
	)

	newton = NewtonIntegrator(damping=0.2)

	O.engines = [
	        ForceResetter(),
	        InsertionSortCollider([Bo1_Sphere_Aabb(), Bo1_Box_Aabb()]),
	        InteractionLoop(
	                [Ig2_Sphere_Sphere_ScGeom(), Ig2_Box_Sphere_ScGeom()], [Ip2_FrictMat_FrictMat_FrictPhys()], [Law2_ScGeom_FrictPhys_CundallStrack()],
	                label="iloop"
	        ),
	        TwoPhaseFlowEngine(dead=1, label="flow"),  #introduced as a dead engine for the moment, see 2nd section
	        GlobalStiffnessTimeStepper(active=1, timeStepUpdateInterval=100, timestepSafetyCoefficient=0.8),
	        triax,
	        newton
	]

	triax.goal1 = triax.goal2 = triax.goal3 = -10000

	while 1:
		O.run(1000, True)
		unb = unbalancedForce()
		if unb < 0.001 and abs(-10000 - triax.meanStress) / 10000 < 0.001:
			break

	setContactFriction(radians(finalFricDegree))

	radius = 0
	for b in O.bodies:
		if b.state.mass == 0:
			b.state.blockedDOFs = 'xyzXYZ'
			b.state.vel = (0, 0, 0)
			b.state.angVel = (0, 0, 0)
		if b.state.mass > 0:
			radius += b.shape.radius
			#b.state.blockedDOFs='xyz'
			#b.state.vel=(0,0,0)
	radius = radius / num_spheres

	triax.dead = True
	while 1:
		O.run(1000, True)
		unb = unbalancedForce()
		if unb < 0.001:
			break

	press = 1000.
	O.run(10, 1)
	flow.dead = 0
	flow.meshUpdateInterval = -1
	flow.useSolver = 3
	flow.permeabilityFactor = 1
	flow.viscosity = 0.1

	flow.bndCondIsWaterReservoir = [0, 0, 1, 0, 0, 0]

	flow.bndCondIsPressure = [0, 0, 1, 0, 0, 0]
	flow.bndCondValue = [0, 0, press, 0, 0, 0]
	flow.boundaryUseMaxMin = [0, 0, 0, 0, 0, 0]
	flow.iniVoidVolumes = True
	newton.damping = 0.1
	GlobalStiffnessTimeStepper.dead = True
	O.dt = builtins.min(0.8 * PWaveTimeStep(), 0.8 * 1. / 1200. * pi / flow.viscosity * graindensity * radius**2)
	O.dynDt = False

	flow.surfaceTension = 0.0
	flow.drainageFirst = False
	flow.isDrainageActivated = False
	flow.isImbibitionActivated = True
	flow.isCellLabelActivated = True
	flow.initialization()
	cs = flow.getClusters()
	c0 = cs[1]

	voidvol = 0.0
	voidvoltot = 0.0
	nvoids = flow.nCells()
	initialvol = [0] * (nvoids)
	bar = [0] * (nvoids)
	initiallevel = O.bodies[flow.wallIds[flow.ymin]
	                       ].state.pos[1] + (O.bodies[flow.wallIds[flow.ymax]].state.pos[1] - O.bodies[flow.wallIds[flow.ymin]].state.pos[1]) / 3

	for ii in range(nvoids):
		initialvol[ii] = 1. / flow.getCellInvVoidVolume(ii)
		voidvoltot += initialvol[ii]
		bar[ii] = flow.getCellBarycenter(ii)[1]

	iniok = 0
	while (iniok == 0):
		celleini1 = [nvoids + 1] * (nvoids)
		celleini0 = [0] * (nvoids)
		for ii in range(len(c0.getInterfaces())):
			if bar[c0.getInterfaces()[ii][1]] < initiallevel:
				if celleini1[c0.getInterfaces()[ii][1]] == nvoids + 1:
					celleini1[c0.getInterfaces()[ii][1]] = ii
					celleini0[c0.getInterfaces()[ii][1]] = c0.getInterfaces()[ii][0]
		for ii in range(nvoids):
			if celleini1[ii] != nvoids + 1:
				flow.clusterOutvadePore(celleini0[ii], ii)
		no = 0
		for ii in range(nvoids):
			if bar[ii] < initiallevel:
				if flow.getCellLabel(ii) == 0:
					no = 1
		if no == 0:
			iniok = 1
			for ii in range(len(c0.getInterfaces())):
				c0.setCapVol(ii, 0.0)

	c0.solvePressure()
	flow.computeCapillaryForce()
	for b in O.bodies:
		O.forces.setPermF(b.id, flow.fluidForce(b.id))
	O.run(1, 1)
	flow.savePhaseVtk("./vtk", True)

	timeini = O.time
	ini = O.iter

	Qin = 0.0
	#Qout=0.0

	totalflux = [0] * (nvoids)
	#totalCellSat=0.0

	for ii in range(nvoids):
		if flow.getCellLabel(ii) == 0:
			voidvol += initialvol[ii]

	bubble = 0
	dd = 0.0
	celleok = [0] * (nvoids)
	deltabubble = 0
	col0 = [0] * (nvoids)
	neighK = [0.0
	         ] * (nvoids)  #FIXME: after remeshing the size will be invalid since nvoids can change, initializations will have to go in the function itself

	def pressureImbibition():
		global Qin, total2, dd, deltabubble, bubble

		c0.updateCapVolList(O.dt)

		Qin += -1. * (flow.getBoundaryFlux(flow.wallIds[flow.ymin])) * O.dt
		#Qout+=(flow.getBoundaryFlux(flow.wallIds[flow.ymax]))*O.dt

		col1 = [0] * (nvoids)
		delta = [0.0] * (nvoids)
		for ii in range(nvoids):
			if flow.getCellLabel(ii) == 0:
				totalflux[ii] += -1. * flow.getCellFluxFromId(ii) * O.dt
				if (totalflux[ii]) >= initialvol[ii]:
					col1[ii] = 1
				if (totalflux[ii]) > initialvol[ii]:
					delta[ii] = totalflux[ii] - initialvol[ii]
					totalflux[ii] += -1 * delta[ii]
					#dd+=delta[ii]

		# advices:
		# never write 'getInterfaces()' inside a loop, it's expensive, get the list once outside loop
		# get interfaces again only if you know the list could have change (cluster got/lost pores).
		# I'm fixing only the first loop below (old version left commented)
		#
		for ii in range(len(c0.getInterfaces())):
			ll = c0.getInterfaces()[ii][1]
			if col1[ll] == 1:
				if celleok[ll] == 0:
					celleok[ll] = 1
					col0[ll] = c0.getInterfaces()[ii][0]

		for jj in range(nvoids):
			if col1[jj] == 1:
				flow.clusterOutvadePore(col0[jj], jj)
				#totalCellSat+=initialvol[jj]

		for ii in range(len(c0.getInterfaces())):
			ll = c0.getInterfaces()[ii][0]
			if delta[ll] != 0:
				neighK[ll] += c0.getConductivity(ii)
		for ii in range(len(c0.getInterfaces())):
			ll = c0.getInterfaces()[ii][0]
			if delta[ll] != 0:
				c0.setCapVol(ii, delta[ll] / neighK[ll] * c0.getConductivity(ii))
				totalflux[c0.getInterfaces()[ii][1]] += delta[ll] / neighK[ll] * c0.getConductivity(ii)
		for ii in range(nvoids):
			if delta[ii] != 0:
				if neighK[ii] == 0:
					deltabubble += delta[ii]
					bubble += 1

		col1 = [0] * (nvoids)
		delta = [0.0] * (nvoids)
		for ii in range(nvoids):
			if flow.getCellLabel(ii) == 0:
				if (totalflux[ii]) >= initialvol[ii]:
					col1[ii] = 1
				if (totalflux[ii]) > initialvol[ii]:
					delta[ii] = totalflux[ii] - initialvol[ii]
					totalflux[ii] += -1 * delta[ii]
					#dd+=delta[ii]
		if col1 != [0] * (nvoids):
			for ii in range(len(c0.getInterfaces())):
				ll = c0.getInterfaces()[ii][1]
				if col1[ll] == 1:
					if celleok[ll] == 0:
						celleok[ll] = 1
						col0[ll] = c0.getInterfaces()[ii][0]
			for jj in range(nvoids):
				if col1[jj] == 1:
					flow.clusterOutvadePore(col0[jj], jj)
					#totalCellSat+=initialvol[jj]
			for ii in range(len(c0.getInterfaces())):
				ll = c0.getInterfaces()[ii][0]
				if delta[ll] != 0:
					neighK[ll] += c0.getConductivity(ii)
			for ii in range(len(c0.getInterfaces())):
				ll = c0.getInterfaces()[ii][0]
				if delta[ll] != 0:
					c0.setCapVol(ii, delta[ll] / neighK[ll] * c0.getConductivity(ii))
					totalflux[c0.getInterfaces()[ii][1]] += delta[ll] / neighK[ll] * c0.getConductivity(ii)
			for ii in range(nvoids):
				if delta[ii] != 0:
					if neighK[ii] == 0:
						deltabubble += delta[ii]
						bubble += 1
			col1 = [0] * (nvoids)
			delta = [0.0] * (nvoids)
			for ii in range(nvoids):
				if flow.getCellLabel(ii) == 0:
					if (totalflux[ii]) >= initialvol[ii]:
						col1[ii] = 1
					if (totalflux[ii]) > initialvol[ii]:
						delta[ii] = totalflux[ii] - initialvol[ii]
						totalflux[ii] += -1 * delta[ii]
						dd += delta[ii]
						print(O.iter, 'waterloss', ii, delta[ii])
			if col1 != [0] * (nvoids):
				for ii in range(len(c0.getInterfaces())):
					ll = c0.getInterfaces()[ii][1]
					if col1[ll] == 1:
						if celleok[ll] == 0:
							celleok[ll] = 1
							col0[ll] = c0.getInterfaces()[ii][0]
				for jj in range(nvoids):
					if col1[jj] == 1:
						flow.clusterOutvadePore(col0[jj], jj)
						#totalCellSat+=initialvol[jj]

		total2 = 0.0
		for ii in range(nvoids):
			total2 += totalflux[ii]

		c0.solvePressure()
		flow.computeCapillaryForce()
		for b in O.bodies:
			O.forces.setPermF(b.id, flow.fluidForce(b.id))

	file = open('Test.txt', "w")
	checkdifference = 0

	def equilibriumtest():
		global F33, F22, checkdifference
		#unbalanced=utils.unbalancedForce()
		F33 = abs(O.forces.f(flow.wallIds[flow.ymax])[1])
		F22 = abs(O.forces.f(flow.wallIds[flow.ymin])[1])
		#F11 =abs(O.forces.f(flow.wallIds[flow.xmax])[0]),
		#F00=abs(O.forces.f(flow.wallIds[flow.xmin])[0]),
		#F44=abs(O.forces.f(flow.wallIds[flow.zmin])[2]),
		#F55=abs(O.forces.f(flow.wallIds[flow.zmax])[2]),
		deltaF = abs(F33 - F22)
		file.write(str(O.iter) + " " + str(F33) + " " + str(F22) + " " + str(deltaF) + "\n")
		if O.time >= timeini + 1.5:
			if checkdifference == 0:
				print('check F done')
				if deltaF > 0.01 * press:
					raise YadeCheckError('Error: too high difference between forces acting at the bottom and upper walls')
					#O.pause()
				checkdifference = 1

	once = 0

	def fluxtest():
		global once, QinOk
		no = 0

		QinOk = Qin - deltabubble
		error = QinOk - total2
		if error > toleranceWarning:
			print(
			        "Warning: difference between total water volume flowing through bottom wall and water loss due to air bubble generations",
			        QinOk, " vs. total water volume flowing inside dry or partially saturated cells", total2
			)
		if error > toleranceCritical:
			raise YadeCheckError("The difference is more, than the critical tolerance!")
		file.write(str(O.time - timeini) + " " + str(total2) + " " + str(QinOk) + " " + str(error) + "\n")

		for ii in range(nvoids):
			if flow.getCellLabel(ii) == 0:
				no = 1
		if once == 0:
			if no == 0:
				imbtime = O.time - timeini
				print(imbtime, voidvol, total2, QinOk)
				if voidvol - total2 > toleranceWarning:
					print(
					        "Warning: initial volume of dry voids", voidvol,
					        " vs. total water volume flowing inside dry or partially saturated cells", total2
					)
				if voidvol - total2 > toleranceCritical:
					raise YadeCheckError("The difference is more, than the critical tolerance!")
				file.write(str(imbtime) + " " + str(voidvol) + " " + str(total2) + " " + str(QinOk) + "\n")
				once = 1
				timing.stats()

	def addPlotData():
		plot.addData(i1=O.iter, t=O.time, Fupper=F33, Fbottom=F22, Q=QinOk, T=total2)

	plot.live = True
	plot.plots = {' t ': ('Fupper', 'Fbottom'), 't': ('Q', 'T')}
	plot.plot()

	def pl():
		flow.savePhaseVtk("./vtk", True)

	O.engines = O.engines + [PyRunner(iterPeriod=100, command='pl()')]
	#O.engines=O.engines+[VTKRecorder(iterPeriod=100,recorders=['spheres'],fileName='./exp')]
	O.engines = O.engines + [PyRunner(iterPeriod=1, command='equilibriumtest()')]
	O.engines = O.engines + [PyRunner(iterPeriod=1, command='pressureImbibition()')]
	O.engines = O.engines + [PyRunner(iterPeriod=1, command='fluxtest()')]
	O.engines = O.engines + [PyRunner(iterPeriod=1, command='addPlotData()')]

	O.timingEnabled = True
	#file.close()
	#plot.saveDataTxt('plots.txt',vars=('i1','t','Fupper','Fbottom','Q','T'))

	#O.run(1,1)

	import tempfile, shutil
	dirpath = tempfile.mkdtemp()
	for fileName in ['./vtk', './Test.txt']:
		if (os.path.exists(fileName)):
			shutil.move(fileName, dirpath)
		print("File %s moved into %s/ directory" % (fileName, dirpath))

else:
	print("This check_TwoPhaseFlowEngine_PressureInjection.py cannot be executed because TWOPHASEFLOW or PFVFLOW are disabled")
