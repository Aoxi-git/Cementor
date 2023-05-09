// 2023 © Karol Brzeziński <karol.brze@gmail.com>
#include <lib/high-precision/Constants.hpp>
#include "SimpleHeatExchanger.hpp"
#include <core/Scene.hpp>

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
	        
	        // as a test put in the test vector 'bId +contactArea(r1 = bId, r2 = bId, penetratioDepth = 0.1*bId)'
	        Real someDummyValue;
	        Real r1, r2, penetrationDepth;
	        r1 = r2 = (Real)bId;
	        penetrationDepth = 0.1*r1;
	        
	        someDummyValue = r1+contactArea(r1,r2,penetrationDepth);
	        test[position] = someDummyValue;
	};	


	return;

}

void SimpleHeatExchanger::action()//
{
    if (previousNumberOfBodies != bodyIds.size()) needsInit = true;
	if (needsInit) init();


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



} // namespace yade
