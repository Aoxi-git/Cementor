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
	
	/*//TESTING CLUMP ID MAPPING
	for (auto clumpCounter = clumpIdtoPosition.cbegin(); clumpCounter != clumpIdtoPosition.cend(); clumpCounter++)// iterating over map
	{
	        vector<long> positions;
	        Body::id_t cId;
	        
	        cId = clumpCounter->first;// key
	        positions = clumpCounter->second;// value
	        
	        long posSize = positions.size();
	        
	    	for (long vCounter = 0; vCounter < posSize; vCounter++) 
	        {
	            long pos = positions[vCounter];
	            test[pos] = (Real)cId;// for test purposes set test value as previously stored cId
	        }
	};*/
	
	//TESTING BODY ID MAPPING + some other functions
	for (auto bodyCounter = bodyIdtoPosition.cbegin(); bodyCounter != bodyIdtoPosition.cend(); bodyCounter++)// iterating over map
	{
	        Body::id_t bId;
	        long position;
	        
	        bId = bodyCounter->first;// key
	        position = bodyCounter->second;// value
	        

	        test[position] = (Real)bId;
	};	


	return;

}

void SimpleHeatExchanger::action()//
{
    if (previousNumberOfBodies != bodyIds.size()) needsInit = true;
	if (needsInit) init();
	
	dTime = scene->time-lastTime;
	lastTime = scene->time;
    energyFlow();
    updateTemp();
    
	previousNumberOfBodies = bodyIds.size();
	return;

}



Real SimpleHeatExchanger::contactArea(Real r1, Real r2, Real penetrationDepth)//Provide radii of both spheres. If one of the radii is 0.0, assume that sphere is contacting facet.
{
	assert(penetrationDepth >= 0);// I don't think this assert works at all.
	assert(r1 >= 0);
	assert(r2 >= 0);
	
	Real d = r1 + r2 - penetrationDepth;//# distance between bodies
	Real a;//radius of intersection circle
	
	if (r1>0 and r2>0)//#two spheres case: https://mathworld.wolfram.com/Sphere-SphereIntersection.html
        a = (0.5/d)*pow((4*pow(d,2)*pow(r1,2)-pow((pow(d,2)-pow(r2,2)+pow(r1,2)),2)),0.5);
    else
        a = pow((2*d*penetrationDepth-pow(penetrationDepth,2)),0.5);//#https://en.wikipedia.org/wiki/Spherical_cap
  
    Real area = Mathr::PI * pow(a,2);
	return area;

}

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
        //#pragma omp parallel for
        for(int j=0; j<nIntr; j++){
           const shared_ptr<Interaction>& i=(*scene->interactions)[j];
           if(!i->isReal()) continue;
		    Body::id_t id1 = i->getId1();
		    Body::id_t id2 = i->getId2();
		    const auto b1  = Body::byId(id1, scene);
		    const auto b2  = Body::byId(id2, scene);
		    Real r1 = 0;// if body is not sphere assume radius = 0
		    Real r2 = 0;
		    
		    shared_ptr<Sphere> sh1 = YADE_PTR_DYN_CAST<Sphere>(b1->shape);//copied from SPherePack.cpp
		    shared_ptr<Sphere> sh2 = YADE_PTR_DYN_CAST<Sphere>(b2->shape);//copied from SPherePack.cpp
		    //const shared_ptr<Sphere>* sh2 = b2->shape.get();
		    if (typeid(sh1) == typeid(Sphere)) r1 = sh1->radius;
		    if (typeid(sh2) == typeid(Sphere)) r2 = sh2->radius;
		    
		    ScGeom*      geom  = dynamic_cast<ScGeom*>(i->geom.get()); 
		    if (typeid(*geom) != typeid(ScGeom)) throw runtime_error("Currently the SimpleHeatExchanger can handle only real interactions of ScGeom type.");
		    Real penetrationDepth = geom->penetrationDepth;
		    
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
    
    Real m1, m2, T1, T2, cond1, cond2, condMin, Eth1, Eth2, EthFlow;
    
    m1 = mass[pos1];
    m2 = mass[pos2];
    T1 = T[pos1];
    T2 = T[pos2];
    cond1 = cond[pos1];
    cond2 = cond[pos2];
    condMin = math::min(cond1, cond2);// conductivity is minimum value of two bodies
    Eth1 = bodyEth[pos1];
    Eth2 = bodyEth[pos2];
    
    //compute how much energy should flow between bodies and check if energy doesn't drop below zero
    EthFlow = condMin*A*dTime*(T1-T2);
    
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
} // namespace yade
