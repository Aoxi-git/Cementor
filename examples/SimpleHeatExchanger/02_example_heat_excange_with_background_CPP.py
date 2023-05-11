# Example based on the gravity deposition example from documentation
# This example, compared to "01_example_CPP.py" contains additional boundary condition 'background'. Uppermost particles (preselected in a simplified way) are cooled down by additional 'dummy' body.
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

#### (NEW PART) now add 'background' for example ambient conditions of given temperature, tak can cool down uppermost particles.
backgroundId = -2 # negative iD for purpose
bodyReal = False
hfcpp.addBody(backgroundId, -1, 0, 0, Tmin, capacity, conductivity, bodyReal)


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



