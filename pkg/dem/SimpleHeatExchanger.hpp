// 2023 © Karol Brzeziński <karol.brze@gmail.com>
#pragma once
#include <core/Scene.hpp>
#include <pkg/common/PeriodicEngines.hpp>


namespace yade { // Cannot have #include directive inside.


class SimpleHeatExchanger : public PeriodicEngine {
private:
	bool  needsInit; // KB: I just leave it temporarily.
	void  init();
	


public:
	
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



//##### BodyTwin

class BodyTwin : public Serializable {
public:
	// numerical types for storing ids
	using id_twin = int;
	// internal structure to hold some interaction of a body; used by InteractionContainer;


	BodyTwin::id_twin   getId() const { return id; };

	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS_CTOR_PY(BodyTwin,Serializable,"Description....",
		((BodyTwin::id_twin,id,,,"Unique id of this body twin."))

		,
		/* ctor */,
		/* py */
		//
	);
	// clang-format on
};
REGISTER_SERIALIZABLE(BodyTwin);

} // namespace yade
