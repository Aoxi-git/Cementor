// 2023 © Karol Brzeziński <karol.brze@gmail.com>
//#pragma once
#include <core/Scene.hpp>
#include <pkg/common/PeriodicEngines.hpp>


namespace yade { // Cannot have #include directive inside.

//##### BodyTwin

class BodyTwin : public Serializable {
public:
	Real Eth;//Thermal energy

    BodyTwin(int id_, Real mass_, Real cap_, Real cond_, Real T_, int clumpId_): id(id_), mass(mass_), T(T_), cap(cap_), cond(cond_), clumpId(clumpId_) {Eth = mass_ * T_ * cap_;};// constructor
    
    Real getEth() const { return Eth;};

	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS_CTOR_PY(BodyTwin,Serializable,"Description....",
		((int,id,,,"Id of this body twin."))
		((Real,mass,,,"Mass of this body twin."))
		((Real,T,,,"T - temperature in [K]."))
		((Real,cap,,,"Specific heat capacity [J/(kg*K)] (449 is value for granite)."))
		((Real,cond,,,"An analog of heat conductivity but the unit is not [W/(m*K)] but [W/(m^2*K)] - need to be found by callibration."))
		((int,clumpId,-1,,"ClumpId of this body twin."))
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
	std::map<int, BodyTwin*> bodyTwins_;
	


public:
	void createBodyTwin(int id_, Real mass_, Real cap_, Real cond_, Real T_, int clumpId_);// nie wiem dlaczego ale ' override' musi być 
	
	void action() override;
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
