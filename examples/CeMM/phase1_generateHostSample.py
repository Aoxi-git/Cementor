# -*- coding: utf-8 -*-
#*************************************************************************
from __future__ import division
from yade import pack, plot
import math
import numpy as np
import random
from random import gauss
import timeit
import pickle

##########  this script is modified from https://gitlab.com/yade-dev/trunk/blob/master/examples/triax-tutorial/script-session1.py#L142
############################################
###   DEFINING VARIABLES AND MATERIALS   ###
############################################
utils.readParamsFromTable(num_spheres=3500)
from yade.params import table

num_spheres=table.num_spheres# number of spheres
targetPorosity = 0.41 #the porosity we want for the packing
compFricDegree = 10 # initial contact friction during the confining phase (will be decreased during the REFD compaction process)
finalFricDegree = 10 # contact friction during the deviatoric loading
damp=0.6 # damping coefficient
stabilityThreshold=0.001 # we test unbalancedForce against this value in different loops (see below)
confinement=100e3

relaxationRatio=1e-8


mn,mx=Vector3(0,0,0),Vector3(0.17,0.17,0.17)  #determine the size of the sample

MatWall=O.materials.append(FrictMat(young=1e10,poisson=0.3,frictionAngle=0,density=0,label='walls'))


MatSand = O.materials.append(CohFrictMat(isCohesive=True,young=1e9,alphaKr=0,alphaKtw=0,\
										 poisson=0.3,frictionAngle=radians(30),etaRoll=0,etaTwist=0,\
										 density=2650.0,normalCohesion=0, shearCohesion=0,\
										 momentRotationLaw=False,label='sand'))


## create walls around the packing
walls=aabbWalls([mn,mx],thickness=0,material='walls')
wallIds=O.bodies.append(walls)

## use a SpherePack object to generate a random loose particles packing
sp=pack.SpherePack()

# =============================================================================
### below can achieve exact psd and num, but not for porosity
# sp.makeCloud(mn,mx,psdSizes=psdSizes,psdCumm=psdCumm,num=num_spheres,distributeMass=True,seed=1)
sp.makeCloud(mn,mx,-1,0.3,num_spheres,False, 0.95,seed=1)
O.bodies.append([sphere(center,rad,material='sand') for center,rad in sp])

Gl1_Sphere.quality=3

newton=NewtonIntegrator(damping=damp)

O.engines=[
	ForceResetter(),
	InsertionSortCollider([Bo1_Sphere_Aabb(),Bo1_Box_Aabb()]),
	InteractionLoop(
		[Ig2_Sphere_Sphere_ScGeom6D(),Ig2_Box_Sphere_ScGeom()],
		[Ip2_CohFrictMat_CohFrictMat_CohFrictPhys(normalCohesion=MatchMaker(matches=((0,0,0),(0,1,0),
																			   (1,1,0))),
												  shearCohesion=MatchMaker(matches=((0,0,0),(0,1,0),
																			   (1,1,0))),
												  label="Ip2Coh"),
		 Ip2_FrictMat_FrictMat_FrictPhys()],
		[Law2_ScGeom6D_CohFrictPhys_CohesionMoment(useIncrementalForm=True,always_use_moment_law=False),Law2_ScGeom_FrictPhys_CundallStrack()]
	),
	GlobalStiffnessTimeStepper(active=1,timeStepUpdateInterval=100,timestepSafetyCoefficient=0.8),
	TriaxialStressController(
		maxMultiplier=1.+3e7/1e9,
		finalMaxMultiplier=1.+16e4/1e9,
		thickness = 0,
		stressMask = 7,
		internalCompaction=False,
		label='triax'),
	TriaxialStateRecorder(iterPeriod=100,file='WallStresses_phase1'),
	newton,
	PyRunner(command='reachConfinement()', iterPeriod=200,label="engineReachConfine"),
	PyRunner(iterPeriod=500,command='reachTargetPorosity()',label='engineTargetPoro'),
]

#Display spheres with 2 colors for seeing rotations better
Gl1_Sphere.stripes=0
#if nRead==0: yade.qt.Controller(), yade.qt.View()


#the value of (isotropic) confining stress defines the target stress to be applied in all three directions
triax.goal1=triax.goal2=triax.goal3=-confinement


engineTargetPoro.dead=True
# engineStopRelaxation.dead=True
def reachConfinement():
	print ('unbF:',unbalancedForce(),' meanStress: ',-triax.meanStress,'top:',-triax.stress(triax.wall_top_id)[1])
	if unbalancedForce() < stabilityThreshold and abs(-confinement-triax.meanStress)/confinement<0.0001:
		O.pause()
		afterReachConfinement()


def afterReachConfinement():
	engineReachConfine.dead=True
	engineTargetPoro.dead=False
	triax.internalCompaction=False
	print ("engineReachConfine closed")
	print ("click run to start reach target porosity")
	O.run()

def reachTargetPorosity():
	global compFricDegree 
	if triax.porosity - targetPorosity>0.0001:
		compFricDegree = 0.95*compFricDegree
		setContactFriction(radians(compFricDegree))
		print ("Friction: ",compFricDegree," porosity:",triax.porosity)
	else:
		O.pause()
		print (" target porosity reached")
		afterReachTargetPorosity()

def afterReachTargetPorosity():
	engineTargetPoro.dead=True
	print ("engineTargetPoro closed")
	setContactFriction(radians(finalFricDegree))
	# NewtonIntegrator.damping=0.9
	# newton.damping=0.9 ## just make sure damping has been changed, according to Yade question:https://answers.launchpad.net/yade/+question/701438

	# engineStopRelaxation.dead=False
	# changeColor()
	# getIniSandIdContactSet()
	# getPickle()
	O.save('HostSample.yade.gz')
	print ("Phase1 completed. Please go to phase2")



# initialSandIds=[] # Using initialSandIds to store ids of sand particles

# def changeColor():
# 	for i in O.bodies:
# 		if isinstance(i.shape,Sphere):
# 			i.shape.color=[1,1,1]
# 			initialSandIds.append(i.id)


# iniSandIdContactSet=set()
# def getIniSandIdContactSet():
# 	for i in O.interactions:
# 		if i.id1 in initialSandIds and i.id2 in initialSandIds:
# 			iniSandIdContactSet.add((i.id1,i.id2))


# # =============================================================================
# def getPickle():
# 	output = open('outputPhase1_100kPa7ksand_young{}_p{}_sandCememtlist.pkl'.format(\
# 				  young,targetPorosity), 'wb')
# 	pickle.dump([initialSandIds,iniSandIdContactSet], output,-1)
# 	# Pickle the list using the highest protocol available.
# 	#pickle.dump(selfref_list, output, -1)
# 	output.close()

utils.waitIfBatch()