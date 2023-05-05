// 2023 © Karol Brzeziński <karol.brze@gmail.com>
//#pragma once
#include <core/Scene.hpp>
#include <pkg/common/PeriodicEngines.hpp>


namespace yade { // Cannot have #include directive inside.

//##### BodyTwin

class BodyTwin : public Serializable {
public:
	Real Eth;//Thermal energy
	
    void  init(int id_t, Real mass_t, Real cap_t, Real cond_t, Real T_t);

	Real  getEth() const { return Eth; };

	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS_CTOR_PY(BodyTwin,Serializable,"Description....",
		((Real,mass,,,"Mass of this body twin."))
		((int,id,,,"Id of this body twin."))
		((int,clumpId,,,"ClumpId of this body twin."))
		((Real,T,,,"T - temperature in [K]."))
		((Real,cap,,,"Specific heat capacity [J/(kg*K)] (449 is value for granite)."))
		((Real,cond,,,"An analog of heat conductivity but the unit is not [W/(m*K)] but [W/(m^2*K)] - need to be found by callibration."))

		,
		/* ctor */,
		/* py */
		//
	);
	// clang-format on
};
REGISTER_SERIALIZABLE(BodyTwin);

class SimpleHeatExchanger : public PeriodicEngine {
private:
	bool  needsInit; // KB: I just leave it temporarily.
	void  init();
	
	


public:
    //typedef std::map<int, BodyTwin*> MyMap;
	using MapId2BodyTwin = std::map<int, shared_ptr<BodyTwin>>;
	
	void action() override;
	//vector<int> twinIds;
	// clang-format off
		YADE_CLASS_BASE_DOC_ATTRS_CTOR(SimpleHeatExchanger,PeriodicEngine,"Description...",
			((vector<int>,twinIds,,,"Ids of BodyTwins."))
			,//!!!!!!!! uwaga - ten przecinek dopiero po wszystkich argumentach
			/*ctor*/ needsInit=true; //  KB: I just leave it temporarily.
		);
	// clang-format on
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(SimpleHeatExchanger);




} // namespace yade
