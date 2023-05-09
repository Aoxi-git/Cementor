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
	return;

}

void SimpleHeatExchanger::action()
{
	if (needsInit) init();
	return;

}





} // namespace yade
