// 2023 © Karol Brzeziński <karol.brze@gmail.com>
#include "SimpleHeatExchanger.hpp"
#include <core/Scene.hpp>

namespace yade { // Cannot have #include directive inside.

YADE_PLUGIN((SimpleHeatExchanger)(BodyTwin));
/************************ UniaxialStrainer **********************/
CREATE_LOGGER(SimpleHeatExchanger);

void BodyTwin::init(int id_t, Real mass_t, Real cap_t, Real cond_t, Real T_t)
{
	id = id_t;
	mass = mass_t;
	cap = cap_t;
	cond = cond_t;
	T = T_t;
	
	Eth = mass * T * cap;
	return;

}

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
