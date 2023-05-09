// 2023 © Karol Brzeziński <karol.brze@gmail.com>
//#pragma once
#include <core/Scene.hpp>
#include <pkg/common/PeriodicEngines.hpp>


namespace yade { // Cannot have #include directive inside.

class SimpleHeatExchanger : public PeriodicEngine {
private:
	bool  needsInit; // KB: I just leave it temporarily.
	void  init();
	std::map<int, int> bodyIdtoPosition;// maps the body id to its position in vector
	


public:
	void action() override;
	// clang-format off
		YADE_CLASS_BASE_DOC_ATTRS_CTOR(SimpleHeatExchanger,PeriodicEngine,"Description...",
		((vector<Real>,mass,,,"Mass of this body twin."))
		((vector<Real>,T,,,"T - temperature in [K]."))
		((vector<Real>,cap,,,"Specific heat capacity [J/(kg*K)] (449 is value for granite)."))
		((vector<Real>,cond,,,"An analog of heat conductivity but the unit is not [W/(m*K)] but [W/(m^2*K)] - need to be found by callibration."))
		((vector<int>,bodyId,,,"Ids of bodies (actual bodies and dummyBodies). It is recommended to use negative values for dummy bodies so it is not mixed with real bodies."))
		((vector<bool>,bodyReal,,,"If true, body is real, else is dummyBody."))
		((vector<Real>,bodyEth,,,"Thermal energy of body. Note it is here for reading purposes only, but I leave it for development phase.."))
		((vector<int>,clumpId,-1,,"ClumpId of this body twin."))
		((vector<Real>,dummyInteractionionA,,,"Areas of interactions."))
		((vector<Real>,dummyInteractionionId1,,,"Id1 of interactions."))
		((vector<Real>,dummyInteractionionId2,,,"Id2 of interactions."))
		((bool,onlyDummyInteractionions,false,,"If true, the heat is exchanged only via dummy interactions."))

			,//!!!!!!!! uwaga - ten przecinek dopiero po wszystkich argumentach
			/*ctor*/ needsInit=true; //  KB: I just leave it temporarily.
		);
	// clang-format on
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(SimpleHeatExchanger);




} // namespace yade
