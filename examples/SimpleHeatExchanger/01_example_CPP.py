# Example based on the gravity deposition example from documentation
from yade import pack

############# BODIES
O.bodies.append(geom.facetBox((.5, .5, .5), (.5, .5, .5), wallMask=31))
sp = pack.SpherePack()
sp.makeCloud((0, 0, 0), (1, 1, 1), rMean=.05, rRelFuzz=.5)
sp.toSimulation()


############# HEAT FLOW SETTINGS
### initialize HeatFlowController
# Note that 'cond' is an analog of head conductivity but the unit is not [W/(m*K)] but [W/(m^2*K)] - need to be found by callibration.
# set temperature of facets temperature to 350 K (except of two side walls)
conductivity = 3e7
capacity = 449.
Tmax = 350
Tmin = 273.15
ids = set([i for i in range(2,10,1)])

hfcpp = SimpleHeatExchanger()
hfcpp.iterPeriod = 100

bodyIds = []
clumpIds = []
mass = []
cap = []
cond = []
T = []
bodyReal = []

for b in O.bodies:
    bodyIds += [b.id]
    clumpIds += [b.clumpId]
    mass += [b.state.mass]
    cap += [capacity]
    cond += [conductivity]
    if b.id in ids:
        T += [Tmax]  
    else:
        T += [Tmin] 
    bodyReal += [True]
    
hfcpp.bodyIds = bodyIds 
hfcpp.clumpIds = clumpIds 
hfcpp.mass = mass
hfcpp.cap = cap
hfcpp.cond = cond 
hfcpp.T = T 
hfcpp.bodyReal = bodyReal  
# colorizing
hfcpp.minT = Tmin
hfcpp.maxT = Tmax  

############ ENGINES
O.engines = [
        ForceResetter(),
        InsertionSortCollider([Bo1_Sphere_Aabb(), Bo1_Facet_Aabb()]),
        InteractionLoop(
                [Ig2_Sphere_Sphere_ScGeom(), Ig2_Facet_Sphere_ScGeom()],
                [Ip2_FrictMat_FrictMat_FrictPhys()],
                [Law2_ScGeom_FrictPhys_CundallStrack()]
        ),
        NewtonIntegrator(gravity=(0, 0, -9.81), damping=0.4),
        hfcpp #PyRunner(command='hf.flow(minT=273,maxT=350)', iterPeriod=100)## heat flow
]
O.dt = .5 * PWaveTimeStep()



