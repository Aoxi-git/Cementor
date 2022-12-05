/*************************************************************************
*  2021 jerome.duriez@inrae.fr                                           *
*  This program is free software, see file LICENSE for details.          *
*************************************************************************/

#ifdef YADE_LS_DEM
#pragma once
#include <core/Dispatching.hpp>
#include <pkg/levelSet/LevelSet.hpp>

namespace yade {
class Bo1_LevelSet_Aabb : public BoundFunctor {
public:
	void go(const shared_ptr<Shape>& cm, shared_ptr<Bound>& bv, const Se3r& se3, const Body*) override;
	FUNCTOR1D(LevelSet);
	YADE_CLASS_BASE_DOC(Bo1_LevelSet_Aabb, BoundFunctor, "Creates/updates an :yref:`Aabb` of a :yref:`LevelSet`");
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(Bo1_LevelSet_Aabb);
} // namespace yade
#endif // YADE_LS_DEM
