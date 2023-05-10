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

#### (NEW PART) now add 'background' for example ambient conditions of given temperature, tak can cool down uppermost particles.
backgroundId = -2 # negative iD for purpose
hfcpp.bodyIds = bodyIds + [backgroundId]
hfcpp.clumpIds = clumpIds + [-1]
hfcpp.mass = mass + [0]# zero mass so the temp would be constant
hfcpp.cap = cap + [capacity]
hfcpp.cond = cond + [conductivity]
hfcpp.T = T + [Tmin]
hfcpp.bodyReal = bodyReal + [False]
# I also need to identify the particles that interactis with this backgroung, I will prepare a function here and will run it separately in Pyrunner.
def create_dummy_interactions():
    """
    I will simply assume that all the particles above z = 0.35 that have less that four interactions should interacti with background
    """
    IDs = []
    for b in O.bodies:
        if isinstance(b.shape,Sphere) and b.state.pos[2]>0.35 and len(b.intrs())<4:
            IDs += [b.id]
    N = len(IDs)
    hfcpp.dummyIntId1 = [backgroundId for i in range(N)]
    hfcpp.dummyIntId2 = IDs
    hfcpp.dummyIntA = [0.00001 for i in range(N)]
    
    
    
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
        PyRunner(command='create_dummy_interactions()', iterPeriod=1000),## update 'uppermost particles'
        hfcpp 
]
O.dt = .5 * PWaveTimeStep()



