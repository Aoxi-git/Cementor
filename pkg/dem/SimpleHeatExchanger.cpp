// 2023 © Karol Brzeziński <karol.brze@gmail.com>
#include "SimpleHeatExchanger.hpp"
#include <core/Scene.hpp>

namespace yade { // Cannot have #include directive inside.

YADE_PLUGIN((SimpleHeatExchanger)(BodyTwin));
/************************ UniaxialStrainer **********************/
CREATE_LOGGER(SimpleHeatExchanger);

/*
void BodyTwin::BodyTwin(int id_, Real mass_, Real cap_, Real cond_, Real T_, int clumpId_)
{
	id = id_;
	mass = mass_;
	cap = cap_;
	cond = cond_;
	T = T_;
	
	Eth = mass * T * cap;
	return;

}*/

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

void SimpleHeatExchanger::createBodyTwin(int id_, Real mass_, Real cap_, Real cond_, Real T_, int clumpId_) 
{
	bodyTwins_[id_] = new BodyTwin(id_, mass_, cap_, cond_, T_, clumpId_);  
	return;

}


} // namespace yade
