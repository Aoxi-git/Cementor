# Example based on the gravity deposition example from documentation
from yade import pack

############# BODIES
O.bodies.append(geom.facetBox((.5, .5, .5), (.5, .5, .5), wallMask=31))
sp = pack.SpherePack()
sp.makeCloud((0, 0, 0), (1, 1, 1), rMean=.05, rRelFuzz=.5)
sp.toSimulation()




############ ENGINES
O.engines = [
        ForceResetter(),
        InsertionSortCollider([Bo1_Sphere_Aabb(), Bo1_Facet_Aabb()]),
        InteractionLoop(
                [Ig2_Sphere_Sphere_ScGeom(), Ig2_Facet_Sphere_ScGeom()],
                [Ip2_FrictMat_FrictMat_FrictPhys()],
                [Law2_ScGeom_FrictPhys_CundallStrack()]
        ),
        NewtonIntegrator(gravity=(0, 0, -9.81), damping=0.4)
]
O.dt = .5 * PWaveTimeStep()


O.run(20000, True)

# after particles drop, divide them in two clumps (on the top) and some spheres touching the bottom. Than initialize the heat flow controller.4
####### divide into clumps

clump1 = []
clump2 = []


for b in O.bodies:
    if isinstance(b.shape,Sphere):
        pos = b.state.pos
        if pos[2] > 0.2: #clump if z> 0.2
            if pos[0] > 0.5:
                clump1 += [b.id]
            else:
                clump2 += [b.id]
            
cId1 = O.bodies.clump(clump1)
cId2 = O.bodies.clump(clump2)       

############# HEAT FLOW SETTINGS
### initialize HeatFlowController
# Note that 'cond' is an analog of head conductivity but the unit is not [W/(m*K)] but [W/(m^2*K)] - need to be found by callibration.
# set temperature of facets temperature to 350 K (except of two side walls)
conductivity = 3e7
capacity = 449.
Tmax = 350
Tmin = 273.15


hfcpp = SimpleHeatExchanger()
hfcpp.iterPeriod = 10
# colorizing
hfcpp.minT = Tmin
hfcpp.maxT = Tmax  

hfcpp.addAllBodiesFromSimulation(Tmin, capacity, conductivity)# first add all bodies
# update properties of some faces
ids = [i for i in range(2,10,1)]
L = [0 for i in range(len(ids))]
hfcpp.addRealBodies(ids, L, Tmax, capacity, conductivity)

### I am not sure what rule could be applied for clumps characteristic length. For now, clumps L need to be updated manually
cIds = [cId1,cId2]
cL = [0.5,0.5]
hfcpp.addRealBodies(cIds, cL , Tmax, capacity, conductivity)

### add engine
O.engines += [hfcpp]

