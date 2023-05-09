// 2023 © Karol Brzeziński <karol.brze@gmail.com>
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

	long        size = bodyIds.size();

	for (long counter = 0; counter < size; counter++) 
	{
		const auto bId = bodyIds[counter];
		const auto cId  = clumpIds[counter];
		test.push_back(0);// prepare proper length of the test vector
		
		bodyIdtoPosition.insert(std::pair<Body::id_t, long>(bId,counter));// map counter to body Ids, so the position in table is stored
		
		if (cId != -1){// store clump info only for clump ids different than -1
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
	}
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
	};

	return;

}

void SimpleHeatExchanger::action()
{
	if (needsInit) init();
	return;

}





} // namespace yade
