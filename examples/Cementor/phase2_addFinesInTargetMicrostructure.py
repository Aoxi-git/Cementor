# -*- coding: utf-8 -*-
#*************************************************************************
from __future__ import division
from yade import pack, plot
import math
import numpy as np
import random
from random import gauss
import timeit



utils.readParamsFromTable(TRmin=0.0003,TRmax=0.00075,\
                          Ncc=15,Tcc=70,\
                          Nco=50,Tco=45,alpha=6.67)
from yade.params import table

TRmin=table.TRmin
TRmax=table.TRmax
Ncc=table.Ncc
Tcc=table.Tcc


Nco=table.Nco
Tco=table.Tco
alpha=table.alpha

O.load('HostSample.yade.gz')
Gl1_Sphere.quality=3


hostSandIds=[] # Using hostSandIds to store ids of sand particles

def getHostSandIds():
    for i in O.bodies:
        if isinstance(i.shape,Sphere):
            hostSandIds.append(i.id)
            i.shape.color=[1,1,1] ## by the way change the color of host sands

getHostSandIds() ## execute now

## define your material for fines, here it is demonstrated as a cohesive material
fineMaterial = O.materials.append(CohFrictMat(isCohesive=True,young=1e9,alphaKr=0,alphaKtw=0,\
                                       poisson=0.3,frictionAngle=radians(19),etaRoll=0,etaTwist=0,\
                                       density=2710.0,normalCohesion=1e6, shearCohesion=1e6,\
                                       momentRotationLaw=True,label='fineMat'))


cohesion=1e10
Ip2Coh.normalCohesion=MatchMaker(matches=((0,0,0),(0,1,0),(0,2,0),(1,1,0),(1,2,cohesion),(2,2,0)))
Ip2Coh.shearCohesion=MatchMaker(matches=((0,0,0),(0,1,0),(0,2,0),(1,1,0),(1,2,cohesion),(2,2,0)))
Ip2Coh.setCohesionNow = True

def fixAllbodies():
    for b in O.bodies:
        if isinstance(b.shape,Sphere):
            b.dynamic=False

fixAllbodies()


def unFixAllbodies():
    for b in O.bodies:
        if isinstance(b.shape,Sphere):
            b.dynamic=True



O.step()


originalParticleContacts=[]

def getOriginalparticleContacts():
    for i in O.interactions:
        if i.id1 in hostSandIds and i.id2 in hostSandIds:
            originalParticleContacts.append(i)

getOriginalparticleContacts() ## execute now


finesIdList_bridging = []
finesIdList_cc = []
finesIdList_co = []

random.seed(40)
selectedOriginalContacts = random.sample(originalParticleContacts, int(len(originalParticleContacts) * Tcc / 100))


selecetedCoatedSandIds=random.sample(hostSandIds,int(len(hostSandIds)*Tco/100))

finesIdSortedByChain = np.zeros((len(selectedOriginalContacts), Ncc), dtype=int)

## Function for generating fine particles in contact cementing type of distribution
def addFines_cc():
    ii=0
    for i in selectedOriginalContacts:
        R1 = O.bodies[i.id1].shape.radius
        R2 = O.bodies[i.id2].shape.radius
        Rmin = min(R1, R2)
        d = i.geom.penetrationDepth ## the overlap of the two particle in contact, delta
        t = R1+R2-d
        if R1 >= R2:
            vectorDirection = -O.bodies[i.id1].state.pos+O.bodies[i.id2].state.pos
            vectorN = vectorDirection/vectorDirection.norm()
            rootargsa = pow(t, 2)-pow(t, 2)/pow(np.sin(np.pi/Ncc),2)-pow(R1-R2,2)
            rootargsb = 2*R1*pow(t, 2)-2*(R1-R2)*(pow(R1,2)+R1*R2-(R1+R2)*d+pow(d,2)/2)
            rootargsc = pow(R1, 2)*pow(t,2)-pow((pow(R1, 2)+R1*R2-(R1+R2)*d+pow(d,2)/2), 2)
            args = [rootargsa, rootargsb, rootargsc]
            fineR = max(np.roots(args))
            h = (R1 - R2) * (fineR + d / 2) / (R1 + R2 - d) 
            ro = sqrt((fineR + R2) ** 2 - (R2 - h - d / 2) ** 2) # ro is the radius of the cement chain.
        else:
            vectorDirection = O.bodies[i.id1].state.pos - O.bodies[i.id2].state.pos
            vectorN = vectorDirection/vectorDirection.norm()
            rootargsa = pow(t, 2) - pow(t, 2) / pow(np.sin(np.pi / Ncc), 2) - pow(R1 - R2, 2)
            rootargsb = 2 * R1 * pow(t, 2) - 2 * (R1 - R2) * (pow(R1, 2) + R1 * R2 - (R1 + R2) * d + pow(d, 2) / 2)
            rootargsc = pow(R1, 2) * pow(t, 2) - pow((pow(R1, 2) + R1 * R2 - (R1 + R2) * d + pow(d, 2) / 2), 2)
            args = [rootargsa, rootargsb, rootargsc]
            fineR = max(np.roots(args))
            h = (R2 - R1) * (fineR + d / 2) / (R1 + R2 - d)
            ro = sqrt((fineR + R1) ** 2 - (R1 - h - d / 2) ** 2)
        # realfineR is a little bigger than fineR, to achieve overlap due to computational truncation
        realfineR = math.ceil((fineR) * 1e5) / 1e5
        a = vectorN[0]
        b = vectorN[1]
        c = vectorN[2]
        xc = i.geom.contactPoint[0]
        yc = i.geom.contactPoint[1]
        zc = i.geom.contactPoint[2]
        cp = np.array([xc,yc,zc]) ## cp is the coordinates of contact point P.
        chainCenter = cp + h*vectorN
        if abs(a) <= abs(b) and abs(a) <= abs(c):
            vectorW=np.array([0, -c, b])
            vectorU = vectorW * ro / np.sqrt(b ** 2 + c ** 2)
            vectorV = np.cross(vectorN, vectorU)
            for k in range(Ncc):
                cf = chainCenter + cos(2*pi*k/Ncc)*vectorU+sin(2*pi*k/Ncc)*vectorV ## cf is the coordinates of a fine particle on that chain.
                finesIdSortedByChain[ii,k]=O.bodies.append(sphere((cf[0], cf[1], cf[2]), radius=realfineR, color=[1, 0, 0],material='fineMat'))
                finesIdList_cc.append(finesIdSortedByChain[ii,k].item())# turn numpy int64 to int
        elif abs(b) <= abs(a) and abs(b) <= abs(c):
            vectorW = np.array([-c, 0, a])
            vectorU = vectorW*ro/np.sqrt(a**2+c**2)
            vectorV = np.cross(vectorN, vectorU)
            for k in range(Ncc):
                cf = chainCenter + cos(2*pi*k/Ncc)*vectorU+sin(2*pi*k/Ncc)*vectorV
                finesIdSortedByChain[ii, k] = O.bodies.append(sphere((cf[0], cf[1], cf[2]), radius=realfineR, color=[1, 0, 0],material='fineMat'))
                finesIdList_cc.append(finesIdSortedByChain[ii, k].item())
        else:
            vectorW = np.array([-b, a, 0])
            vectorU = vectorW * ro / np.sqrt(a ** 2 + b ** 2)
            vectorV = np.cross(vectorN, vectorU)
            for k in range(Ncc):
                cf = chainCenter + cos(2*pi*k/Ncc)*vectorU+sin(2*pi*k/Ncc)*vectorV
                finesIdSortedByChain[ii, k] = O.bodies.append(sphere((cf[0], cf[1], cf[2]), radius=realfineR, color=[1, 0, 0],material='fineMat'))
                finesIdList_cc.append(finesIdSortedByChain[ii, k].item())
        ii+=1
        

## Function for generating fine particles in bridging type of distribution
def addFines_Briging():
    global TRmin,TRmax
    for i in hostSandIds:
        for j in hostSandIds:
            if i>=j:
                continue
            else:
                R1=O.bodies[i].shape.radius
                R2=O.bodies[j].shape.radius
                vectorS1=O.bodies[i].state.pos
                vectorS2=O.bodies[j].state.pos
                centerDist = (vectorS2 - vectorS1).norm()
                d = centerDist - (R1 + R2)
                fineR=d/2.0
                realfineR = math.ceil((d / 2.0) * 1e7) / 1e7
                Rmin = min(R1, R2)
                Rmax = max(R1,R2)
                upLimit = TRmax
                lowLimit = TRmin
                vectorN = (vectorS2 - vectorS1) / centerDist
                GapCenter = (0.5*d + R1) * vectorN + vectorS1

                if realfineR >= upLimit:
                    pass
                elif realfineR < lowLimit:
                    pass
                else:
                    finesIdList_bridging.append(
                        O.bodies.append(sphere((GapCenter[0], GapCenter[1], GapCenter[2]), radius=realfineR, color=[0, 0.8, 0],material='fineMat')))



### function for generating random unit vector
def randomUnitVector():
    vec = [gauss(0, 1) for i in range(3)]
    mag = sum(x**2 for x in vec) ** .5
    return [x/mag for x in vec]

## Function for generating fine particles in coating type of distribution
def addFines_coating():
    global alpha
    for j in selecetedCoatedSandIds:
        sandR = O.bodies[j].shape.radius
        fineR = sandR/alpha
        realfineR = math.ceil(fineR * 1e9) / 1e9
        ## dists is the distance between sand center and fine center
        dists=sandR+fineR
        sandCenter=O.bodies[j].state.pos
        ### generate a random unit vector which is stored in this list, if we don't store it, then, the components in 3 directions are not matched
        randomList = []
        for i in range(Nco):
            randomList.append(randomUnitVector())
            ## the vector from sand center to Ca center
            vectorSand_fine=dists*np.array([randomList[i][0],randomList[i][1],randomList[i][2]])
            fineCenter=sandCenter+vectorSand_fine
            finesIdList_co.append(O.bodies.append(sphere((fineCenter[0], fineCenter[1], fineCenter[2]), radius=realfineR, color=[0, 0, 0.8],material='fineMat')))


undesiredFinesId_bridging=set() ## Store cement particle id

undesiredFinesId_cc=set() 
undesiredFinesId_co=set() 
undesiredFinesId_co=set()
undesiredFinesId_total=set() ## sum id of cements which should be deleted

finalFinesId_bridging=set()
finalFinesId_cc=set()   ## Store id of ultimately left contact cementing particles, after erasing overlap carbonate
finalFinesId_co=set()
finalFinesId_total=set() 


def findUndesiredFines():
    global undesiredFinesId_cc,undesiredFinesId_co,undesiredFinesId_bridging
    ## bridging, coating should not have fine-fine contacts
    for i in O.interactions:
        if i.id1 in finesIdList_bridging and i.id2 in finesIdList_bridging:
            undesiredFinesId_bridging.add(i.id1)
        if i.id1 in finesIdList_co and i.id2 in finesIdList_co:
            undesiredFinesId_co.add(i.id1)
    ## all cases should not have fine-wall contacts
    for i in finesIdList_cc:
        for j in range(6):
            if O.interactions.has(i,j):
                if O.interactions[i, j].isReal:
                    undesiredFinesId_cc.add(i)
    for i in finesIdList_co:
        for j in range(6):
            if O.interactions.has(i,j):
                if O.interactions[i, j].isReal:
                    undesiredFinesId_co.add(i)
    for i in finesIdList_bridging:
        for j in range(6):
            if O.interactions.has(i,j):
                if O.interactions[i, j].isReal:
                    undesiredFinesId_bridging.add(i)
     ## double check according to interaction number 
     ## fines of bridging should have 2 interactions, 
     ## fines of contact cementing is 4, coating is 1
    for i in finesIdList_cc:
        if len(O.bodies[i].intrs()) != 4:
            undesiredFinesId_cc.add(i)
    for i in finesIdList_bridging:
        if len(O.bodies[i].intrs()) != 2:
            undesiredFinesId_bridging.add(i)
    for i in finesIdList_co:
        if len(O.bodies[i].intrs()) != 1:
            undesiredFinesId_co.add(i)


## ## undesiredFinesId_cc is the individual cement particles, and based on undesiredFinesId_cc, we can find the non-closed cement chain which are deleted if you want to make a closed chain.
undesiredFinesId_cc_extendedByChain=set()
def findUndesiredFines_cc_extendedByChain():
    for i in undesiredFinesId_cc:
        loc=np.where(finesIdSortedByChain==i)
        for j in finesIdSortedByChain[loc[0][0]]:
            undesiredFinesId_cc_extendedByChain.add(j.item())


### Delete non-wanted cement particles
def removeUndesiredFines():
    global undesiredFinesId_cc,undesiredFinesId_co,undesiredFinesId_bridging,undesiredFinesId_total
    undesiredFinesId_total=undesiredFinesId_cc_extendedByChain|undesiredFinesId_co|undesiredFinesId_bridging
    for i in undesiredFinesId_total:
        O.bodies.erase(i)

##undesiredFinesId_cc_extendedByChain is int64,
## undesiredFinesId_co, undesiredFinesId_bridging are int

def findFinalFines():
    global finalFinesId_cc,finalFinesId_co,finalFinesId_bridging,finalFinesId_total
    for i in finesIdList_cc:
        finalFinesId_cc.add(i)
    for i in undesiredFinesId_cc_extendedByChain:
        finalFinesId_cc.remove(i)
    for i in finesIdList_co:
        finalFinesId_co.add(i)
    for i in undesiredFinesId_co:
        finalFinesId_co.remove(i)
    for i in finesIdList_bridging:
        finalFinesId_bridging.add(i)
    for i in undesiredFinesId_bridging:
        finalFinesId_bridging.remove(i)
    finalFinesId_total=finalFinesId_cc|finalFinesId_co|finalFinesId_bridging


def mixSample():
    addFines_Briging()
    addFines_cc()
    addFines_coating()
    O.step()
    Ip2Coh.setCohesionNow=True
    findUndesiredFines()
    findUndesiredFines_cc_extendedByChain()
    removeUndesiredFines()
    findFinalFines()
    O.save('TRmin{}_TRmax{}_Tcc{}_Tco{}.yade.gz'.format(TRmin,TRmax,Tcc,Tco))

mixSample()

# =============================================================================

def getEachTypeCementMassContentRatio():
    mass_cc=mass_coat=mass_brig=mass_totalCement=0 # mass of fines in each distribution pattern, cc means contact cementing.
    ratioInTotalCement_cc=ratioInTotalCement_coating=0
    ratioInTotalCement_bridging=0
    global finalFinesId_cc,finalFinesId_co,finalFinesId_bridging,\
    finalFinesId_total
    for i in finalFinesId_cc:
        mass_cc +=O.bodies[i].state.mass
    for i in finalFinesId_co:
        mass_coat +=O.bodies[i].state.mass
    for i in finalFinesId_bridging:
        mass_brig +=O.bodies[i].state.mass
    # for i in CaBeingLeft_poreFill:
    #     mass_poreFill +=O.bodies[i].state.mass
    mass_totalCement=mass_cc+mass_coat+mass_brig
    ratioInTotalCement_cc=mass_cc/mass_totalCement
    ratioInTotalCement_coat=mass_coat/mass_totalCement
    ratioInTotalCement_brig=mass_brig/mass_totalCement
    # ratioInTotalCement_poreFill=mass_poreFill/mass_totalCement
    return 'mass_cc = ',mass_cc,'mass_bridging = ', mass_brig,\
    'mass_coating = ',mass_coat,'mass_totalCement = ',mass_totalCement,\
    'ratioInTotalCement_cc = ',ratioInTotalCement_cc,'ratioInTotalCement_bridging = ',ratioInTotalCement_brig,\
           'ratioInTotalCement_coating = ',ratioInTotalCement_coat


def getCementMassContent():
    totalMassOfSand =totalMassOfCement =cementMassContent =0
    for i in finalFinesId_total:
        totalMassOfCement += O.bodies[i].state.mass
    for i in hostSandIds:
        totalMassOfSand += O.bodies[i].state.mass
    cementMassContent=totalMassOfCement/(totalMassOfSand+totalMassOfCement)
    cs=totalMassOfCement/totalMassOfSand
    return 'totalMassOfSand = ',totalMassOfSand,'totalMassOfCement = ',\
    totalMassOfCement,'cementMassContent = ',cementMassContent,'cementMass/SandMass = ',cs


## A function for calculating coordination number of sands and fines respectively.
def coordinationNumber():
    CN_s = 0  # contact number of sand
    CN_f = 0  # contact number of fine
    for i in hostSandIds:
        for j in O.bodies[i].intrs():
            if j.isReal:
                CN_s += 1
    for j in finalFinesId_total:
        for k in O.bodies[j].intrs():
            if k.isReal:
                CN_f += 1
    CN_all = CN_s + CN_f
    Z_s= CN_s / len(hostSandIds) # coordination number of sand
    Z_f= CN_f / len(finalFinesId_total) # coordination number of fine
    Z_all = CN_all/(len(hostSandIds)+len(finalFinesId_total)) # coordination number of all particles
    return Vector3(Z_s,Z_f,Z_all)

    


## Here you may output some variables if you are interested.
f = open("outputCementedSampleInfor_TRmin{}_TRmax{}_Tcc{}_Tco{}_alpha{}.txt".format(TRmin,TRmax,Tcc,Tco,alpha), "a+")
f.write('Input params: TRmin = {},TRmax = {},Ncc = {},Tcc = {},Nco = {},Tco = {},alpha = {} \n' .format\
            (TRmin,TRmax,Ncc,Tcc,Nco,Tco,alpha))

f.write('{} \n' .format(getCementMassContent()))

f.write('{} \n' .format(getEachTypeCementMassContentRatio()))

f.write('Total number of fines left = {},\n number of fines in contact cementing pattern left = {},\n number of fines in bridging pattern left = {},\n number of fines in coaing pattern left= {} \n' .format\
            (len(finalFinesId_total),len(finalFinesId_cc),len(finalFinesId_bridging),len(finalFinesId_co)))

f.write('Coordination number of sand = {},\n Coordination number of fine = {},\n Coordination number of all particles = {},\n ' .format\
            (coordinationNumber()[0],coordinationNumber()[1],coordinationNumber()[2]))

f.close()


############
utils.waitIfBatch()