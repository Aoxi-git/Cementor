// 2023 © Karol Brzeziński <karol.brze@gmail.com>
#include <lib/high-precision/Constants.hpp>
#include "SimpleHeatExchanger.hpp"
#include "ScGeom.hpp"
#include <core/Scene.hpp>
#include <pkg/common/Sphere.hpp>

namespace yade { // Cannot have #include directive inside.

YADE_PLUGIN((SimpleHeatExchanger));
/************************ UniaxialStrainer **********************/
CREATE_LOGGER(SimpleHeatExchanger);

/*############ INITIALIZATION ##############*/
void SimpleHeatExchanger::init()
{
	needsInit = false;
	assert(bodyIds.size() > 0);
	assert(bodyIds.size() == clumpIds.size());
	
    // clear initialized values
    bodyIdtoPosition.clear();
    clumpIdtoPosition.clear();
    bodyEth.clear();
    test.clear();	

	long        size = bodyIds.size();

	for (long counter = 0; counter < size; counter++) 
	{
		const auto bId = bodyIds[counter];
		const auto cId  = clumpIds[counter];
		test.push_back(0);// prepare proper length of the test vector
		// Prepare map of body IDs
		if (bodyIdtoPosition.count(bId) == 0 ) 
		{
		    bodyIdtoPosition.insert(std::pair<Body::id_t, long>(bId,counter));// map counter to body Ids, so the position in table is stored
		} 
		else
		{
            throw runtime_error("All bodyIds should be unique");
		};
		
		// Prepare map of clump IDs
		if (cId != -1)
		{// store clump info only for clump ids different than -1
		    if (clumpIdtoPosition.count(cId) == 0 ) //if there is no current cId in the map (if count>0 is interpreted as true)
		    {
		        vector<long> positions;
		        positions.push_back(counter);
		        clumpIdtoPosition.insert(std::pair<Body::id_t, vector<long>>(cId,positions));
		    } 
		    else 
		    {
		         vector<long> positions;
		         positions = clumpIdtoPosition[cId];
		         positions.push_back(counter);
		         clumpIdtoPosition[cId] = positions;
		    };
		};
		// Initialize Eth (compute Eth based on Temperature in the initialization phase. Usually compute Temperature based on Eth)
		Real Eth;
		Eth = mass[counter] * T[counter] * cap[counter];
		bodyEth.push_back(Eth);// prepare proper length of the bodyEth vector
	}

	return;

}

/*############ BODIES AND INTERACIONS HANDLING ##############*/

void SimpleHeatExchanger::addRealBody(Body::id_t bId, Real L_, Real T_, Real cap_, Real cond_)  
{
    if (bodyIdtoPosition.size() == 0) init();// initialize if never initialized
    needsInit = true;// also set flag for initialization since the body setup changed
    
    const auto b  = Body::byId(bId, scene);
    if (!b) throw runtime_error("There is no body with given Id in the simulation.");
    Real m = b->state->mass;
    Body::id_t cId = b->clumpId;
    
    addBody(bId, cId, m, L_, T_, cap_, cond_, true);
	
	return;
}

void SimpleHeatExchanger::addBody(Body::id_t bId, Body::id_t cId, Real mass_, Real L_, Real T_, Real cap_, Real cond_, bool Real_)
{
    if (bodyIdtoPosition.size() == 0) init();// initialize if never initialized
    needsInit = true;// also set flag for initialization since the body setup changed   
    // If body id is in bodyIds, update properties, otherwise add Id at the end.
    if (bodyIdtoPosition.count(bId) == 0 )
    {
        bodyIds.push_back(bId);
        mass.push_back(mass_);
        L.push_back(L_);
        T.push_back(T_);
        cap.push_back(cap_);
        cond.push_back(cond_);
        bodyReal.push_back(Real_);
        clumpIds.push_back(cId);        
    }
    else
    {
        long pos = bodyIdtoPosition[bId];
        mass[pos] = mass_;
        L[pos] = L_;
        T[pos] = T_;
        cap[pos] = cap_;
        cond[pos] = cond_;
        bodyReal[pos] = Real_;
        clumpIds[pos] = cId;
    }
	
	return;
}



void SimpleHeatExchanger::addAllBodiesFromSimulation(Real T_ , Real cap_, Real cond_ )
{
    long nBodies = scene->bodies->size();
    
	for (long counter = 0; counter < nBodies; counter++) 
	{
	    const shared_ptr<Body>& b=(*scene->bodies)[counter];
	    Body::id_t bId = b->id;
	    Real L_ = 0;
	    
	    Sphere* sh = dynamic_cast<Sphere*>(Body::byId(bId)->shape.get());// Based on VTKRecorder.cpp and  SPherePack.cpp
	    if (sh) L_ = sh->radius; //
	    
	    addRealBody(bId, L_, T_, cap_, cond_);
	};     
	init();
	return;
}

void SimpleHeatExchanger::addRealBodies(vector<Body::id_t> bodyIds_,  vector<Real> vectL_, Real T_, Real cap_, Real cond_)
{
    long        size = bodyIds_.size();
	for (long counter = 0; counter < size; counter++) 
	{
	    Body::id_t bId = bodyIds_[counter];
	    Real L_ = vectL_[counter];
	    addRealBody(bId, L_, T_, cap_, cond_);
	};  

	init();
	return;
}
/*############ MAIN ACTION ##############*/

void SimpleHeatExchanger::action()//
{
    if (previousNumberOfBodies != bodyIds.size()) needsInit = true;
	if (needsInit) init();
	
	dTime = scene->time-lastTime;
	lastTime = scene->time;
    energyFlow();
    updateTemp();
    
    if (colorize) updateColors();
    
	previousNumberOfBodies = bodyIds.size();
	return;

}


/*############ HEAT FLOW AND TEMPERATURE CONTROLL ##############*/


void SimpleHeatExchanger::energyFlow()//
{
    /* Let energy flow based on the temperature differences between interacting bodies.*/
    // first dummy interactions
    long size = dummyIntId1.size();
    for (long i = 0; i < size; i++)
    {
         Body::id_t id1 = dummyIntId1[i];   
         Body::id_t id2 = dummyIntId2[i];   
         Real A = dummyIntA[i];
         energyFlowOneInteraction(id1, id2, A);
    };
    
    if (!onlyDummyInt)
    {

        int nIntr=(int)scene->interactions->size(); // hoist container size
#ifdef YADE_OPENMP
#pragma omp parallel for schedule(guided) num_threads(ompThreads > 0 ? std::min(ompThreads, omp_get_max_threads()) : omp_get_max_threads())
#endif
        for(int j=0; j<nIntr; j++){
           const shared_ptr<Interaction>& i=(*scene->interactions)[j];
           if(!i->isReal()) continue;
		    Body::id_t id1 = i->getId1();
		    Body::id_t id2 = i->getId2();
		    const auto b1  = Body::byId(id1, scene);
		    const auto b2  = Body::byId(id2, scene);
		    Real r1 = 0;// if body is not sphere assume radius = 0
		    Real r2 = 0;

		    
		    ScGeom*      geom  = dynamic_cast<ScGeom*>(i->geom.get()); 
		    if (typeid(*geom) != typeid(ScGeom)) throw runtime_error("Currently the SimpleHeatExchanger can handle only real interactions of ScGeom type.");
		    Real penetrationDepth = geom->penetrationDepth;

	        Sphere* sh1 = dynamic_cast<Sphere*>(Body::byId(id1)->shape.get());// Based on VTKRecorder.cpp and  SPherePack.cpp
	        if (sh1) r1 = sh1->radius; //
	        Sphere* sh2 = dynamic_cast<Sphere*>(Body::byId(id2)->shape.get());// Based on VTKRecorder.cpp and  SPherePack.cpp
	        if (sh2) r1 = sh2->radius; //

		    
		    Real A = contactArea(r1, r2, penetrationDepth);
		    energyFlowOneInteraction(id1, id2, A);
        }        


    }
    
	return;

}

void SimpleHeatExchanger::energyFlowOneInteraction(Body::id_t id1, Body::id_t id2, Real A)//
{
    /* Let energy flow based on the temperature differences between interacting bodies. Version for one interaction to optimize the energyFlow() function.*/
    if (A < 0 ) throw runtime_error("Area cannot be negative.");
    
    long pos1 = bodyIdtoPosition[id1]; //position of body1 data in all the vectors
    long pos2 = bodyIdtoPosition[id2];
    
    Body::id_t cId1 = clumpIds[pos1];
    Body::id_t cId2 = clumpIds[pos2];
    
    
    if (cId1 != -1) pos1 = bodyIdtoPosition[cId1]; // If body is clumped, threat the whole clump as the the energy source
    if (cId2 != -1) pos2 = bodyIdtoPosition[cId2]; 
    
    Real m1, m2, L1, L2, T1, T2, cond1, cond2, condMean, Eth1, Eth2, EthFlow;
    
    m1 = mass[pos1];
    m2 = mass[pos2];
    L1 = L[pos1];
    L2 = L[pos2];
    if (L1 == 0 and L2 == 0) throw runtime_error("Characteristic lengtch (L) of two bodies exchanging heat cannot be equal to zero.");
    T1 = T[pos1];
    T2 = T[pos2];
    cond1 = cond[pos1];
    cond2 = cond[pos2];
    condMean = 1/(L1/cond1+L2/cond2);// heat pipe concept
    Eth1 = bodyEth[pos1];
    Eth2 = bodyEth[pos2];
    
    //compute how much energy should flow between bodies and check if energy doesn't drop below zero
    EthFlow = condMean*A*dTime*(T1-T2);
    
    if (Eth1 - EthFlow < 0 and m1 > 0) //I check mass because zero mass mean the body temperature is constant. So if m1 = 0 its energy is unlimited.
    {
        EthFlow = Eth1;// maximum energy flow is limited to the energy stored in body 1
    }
    else if (Eth2 + EthFlow < 0 and m2 > 0)
    {
        EthFlow = -Eth2; // minus takes into account that must be in case of T2>T1 and consequently flow from body 2 to body 1
    }

    // modify energy of the body, only if mass > 0
    if (m1 > 0)
        bodyEth[pos1] = Eth1-EthFlow;
    if (m2 > 0)
        bodyEth[pos2] = Eth2+EthFlow;
	return;

}



void SimpleHeatExchanger::updateTemp()//
{
    /* Update temperature of the bodies based on their energy. Match the temperature of the clump members to the temperature of the clump.*/
    long size = bodyIds.size();
#ifdef YADE_OPENMP
#pragma omp parallel for schedule(guided) num_threads(ompThreads > 0 ? std::min(ompThreads, omp_get_max_threads()) : omp_get_max_threads())
#endif    
    for (long i = 0; i < size; i++)
    {
        Real m = mass[i];// Here I don't need position mapping since I iterate over bodyIds
        Real Eth = bodyEth[i];
        Real bCap = cap[i];
        Body::id_t bId = bodyIds[i];
        Body::id_t cId = clumpIds[i];
        Real Temp = Eth/(m * bCap);// 
        
        if (m > 0) // Only change temperature of bodies with mass
        {
            if (cId == -1)// If Standalone
            {
                T[i] = Temp;
            } 
            else if (bId == cId) // is Clump
            {
                T[i] = Temp;// update clump T
                // update T of clump members
                vector<long> positions;// positions of the clump members
                positions = clumpIdtoPosition[cId];
                long posSize = positions.size();
	        
	        	for (long vCounter = 0; vCounter < posSize; vCounter++) 
	            {
	                long pos = positions[vCounter];
	                T[pos] = Temp;
	            }
                
            }// Else is not necessary since clump members were handled with the clump
            
        }
    };  
}

/*############ PRE AND POST PROCESSING ##############*/
Real SimpleHeatExchanger::contactArea(Real r1, Real r2, Real penetrationDepth)//Provide radii of both spheres. If one of the radii is 0.0, assume that sphere is contacting facet.
{
	assert(penetrationDepth >= 0);// I don't think this assert works at all.
	assert(r1 >= 0);
	assert(r2 >= 0);
	if (r1 == 0 and r2 == 0) throw runtime_error("Both radii cannot be equal to zero for ScGeom.");
	
	Real d = r1 + r2 - penetrationDepth;//# distance between bodies
	Real a;//radius of intersection circle
	
	if (r1>0 and r2>0)//#two spheres case: https://mathworld.wolfram.com/Sphere-SphereIntersection.html
        a = (0.5/d)*pow((4*pow(d,2)*pow(r1,2)-pow((pow(d,2)-pow(r2,2)+pow(r1,2)),2)),0.5);
    else
        a = pow((2*d*penetrationDepth-pow(penetrationDepth,2)),0.5);//#https://en.wikipedia.org/wiki/Spherical_cap
  
    Real area = Mathr::PI * pow(a,2);
	return area;

}

void SimpleHeatExchanger::updateColors()//
{
    /* Update temperature of the bodies based on their energy. Match the temperature of the clump members to the temperature of the clump.*/
    long size = bodyIds.size();
#ifdef YADE_OPENMP
#pragma omp parallel for schedule(guided) num_threads(ompThreads > 0 ? std::min(ompThreads, omp_get_max_threads()) : omp_get_max_threads())
#endif    
    for (long i = 0; i < size; i++)
    {
        bool bReal = bodyReal[i];// Here I don't need position mapping since I iterate over bodyIds
        if (bReal)
        {
            Real Temp = T[i];
            Real red, green, blue;
            //red
            red = (Temp-minT)/(maxT-minT);
            //blue
            if((Temp-minT)/(maxT-minT)<0.5)
            {
                blue = 1-2*(Temp-minT)/(maxT-minT) ;
            }
            else
            {
                blue = 0;
            }
            //green
            if((Temp-minT)/(maxT-minT)<0.75)
            {
                green = 0;
            }
            else
            {
                green = 4*(Temp-minT)/(maxT-minT)-3;
            }
                
            Body::id_t id = bodyIds[i];
            auto b  = Body::byId(id, scene);
		    if (!b) continue;
		    Vector3r color = Vector3r(red,green,blue);
		    b->shape->color = color;
		    
            
		}
    }; 
    return; 
}

} // namespace yade
