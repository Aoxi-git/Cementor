
#pragma once
#include <core/Dispatching.hpp>
#include <core/IGeom.hpp>
#include <core/IPhys.hpp>
#include <core/Shape.hpp>
#include <core/State.hpp>
#include <pkg/common/Facet.hpp>
#include <pkg/common/Sphere.hpp>
#include <pkg/common/Wall.hpp>
#include <pkg/dem/DemXDofGeom.hpp>
#include <pkg/dem/FrictPhys.hpp>
#ifdef YADE_OPENGL
#include <pkg/common/GLDrawFunctors.hpp>
#endif

namespace yade { // Cannot have #include directive inside.

/*

L3Geom Ig2 functor cooperation:

1. functors define genericGo function, which take ``bool is6Dof`` as the first arg; they are called from functors returning both L3Geom and L6Geom (this only applies to sphere+sphere now, since Ig2_{Facet,Wall}_Sphere_L6Geom has not been written yet for lack of interest, though would be trivial)
2. genericGo function computes several parameter specific to shape types; then it calls (if L3GEOM_SPHERESLIKE is defined) handleSpheresLikeContact which contains the common code, given all data necessary.

Note that:

(a) although L3GEOM_SPHERESLIKE is enabled by default, its performance impact has not been measured yet (the compiler should be smart enough, since it is just factoring out common code).
(b) L3Geom only contains contPt and normal, which supposes (in L3Geom::applyLocalForce) that particles' centroids and the contact point are colinear; while this is true for spheres and mostly OK for facets&walls (since they are non-dynamic), it might be adjusted in the future -- L3Geom_Something deriving from L3Geom will be created, and exact branch vectors contained in it will be gotten via virtual method from L3Geom::applyLocalForce. This would be controlled via some approxMask (in the Law2 functor perhaps) so that it is only optional
(c) Ig2_Facet_Sphere_L3Geom is only enabled with L3GEOM_SPHERESLIKE

*/

#define L3GEOM_SPHERESLIKE

struct L3Geom : public GenericSpheresContact {
	virtual ~L3Geom();

	// utility function
	// TODO: currently supposes body's centroids are conencted with distance*normal
	// that will not be true for sphere+facet and others, watch out!
	// the force is oriented as applied to particle #1
	void applyLocalForce(const Vector3r& f, const Interaction* I, Scene* scene, NormShearPhys* nsp = NULL) const;
	void applyLocalForceTorque(const Vector3r& f, const Vector3r& t, const Interaction* I, Scene* scene, NormShearPhys* nsp = NULL) const;

	Vector3r relU() const { return u - u0; }

	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS_INIT_CTOR_PY(L3Geom,GenericSpheresContact,"Geometry of contact given in local coordinates with 3 degress of freedom: normal and two in shear plane. [experimental]",
		((Vector3r,u,Vector3r::Zero(),,"Displacement components, in local coordinates. |yupdate|"))
		((Vector3r,u0,Vector3r::Zero(),,"Zero displacement value; u0 should be always subtracted from the *geometrical* displacement *u* computed by appropriate :yref:`IGeomFunctor`, resulting in *u*. This value can be changed for instance\n\n#. by :yref:`IGeomFunctor`, e.g. to take in account large shear displacement value unrepresentable by underlying geomeric algorithm based on quaternions)\n#. by :yref:`LawFunctor`, to account for normal equilibrium position different from zero geometric overlap (set once, just after the interaction is created)\n#. by :yref:`LawFunctor` to account for plastic slip.\n\n.. note:: Never set an absolute value of *u0*, only increment, since both :yref:`IGeomFunctor` and :yref:`LawFunctor` use it. If you need to keep track of plastic deformation, store it in :yref:`IPhys` isntead (this might be changed: have *u0* for :yref:`LawFunctor` exclusively, and a separate value stored (when that is needed) inside classes deriving from :yref:`L3Geom`."))
		/* Is it better to store trsf as Matrix3 or Quaternion?
		* Quaternions are much easier to re-normalize, which we should do to avoid numerical drift.
		* Multiplication of vector with quaternion is internally done by converting to matrix first, anyway
		* We need to extract local axes, and that is easier to be done from Matrix3r (columns)
		*/
		((Matrix3r,trsf,Matrix3r::Identity(),,"Transformation (rotation) from global to local coordinates. (the translation part is in :yref:`GenericSpheresContact.contactPoint`)"))
		((Vector3r,F,Vector3r::Zero(),,"Applied force in local coordinates [debugging only, will be removed]"))
		,
		/*init*/
		,
		/*ctor*/ createIndex();
		, /*py*/
	);
	// clang-format on
	REGISTER_CLASS_INDEX(L3Geom, GenericSpheresContact);
};
REGISTER_SERIALIZABLE(L3Geom);

struct L6Geom : public L3Geom {
	virtual ~L6Geom();
	Vector3r relPhi() const { return phi - phi0; }
	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS_CTOR(L6Geom,L3Geom,"Geometric of contact in local coordinates with 6 degrees of freedom. [experimental]",
		((Vector3r,phi,Vector3r::Zero(),,"Rotation components, in local coordinates. |yupdate|"))
		((Vector3r,phi0,Vector3r::Zero(),,"Zero rotation, should be always subtracted from *phi* to get the value. See :yref:`L3Geom.u0`."))
		,
		/* ctor */ createIndex();
	);
	// clang-format on
	REGISTER_CLASS_INDEX(L6Geom, L3Geom);
};
REGISTER_SERIALIZABLE(L6Geom);

#ifdef YADE_OPENGL
struct Gl1_L3Geom : public GlIGeomFunctor {
	RENDERS(L3Geom);
	void go(const shared_ptr<IGeom>&, const shared_ptr<Interaction>&, const shared_ptr<Body>&, const shared_ptr<Body>&, bool) override;
	void draw(const shared_ptr<IGeom>&, bool isL6Geom = false, const Real& phiScale = 0);
	// clang-format off
	YADE_CLASS_BASE_DOC_STATICATTRS(Gl1_L3Geom,GlIGeomFunctor,"Render :yref:`L3Geom` geometry.",
		((bool,axesLabels,false,,"Whether to display labels for local axes (x,y,z)"))
		((Real,axesScale,1.,,"Scale local axes, their reference length being half of the minimum radius."))
		((Real,axesWd,1.,,"Width of axes lines, in pixels; not drawn if non-positive"))
		((Real,uPhiWd,2.,,"Width of lines for drawing displacements (and rotations for :yref:`L6Geom`); not drawn if non-positive."))
		((Real,uScale,1.,,"Scale local displacements (:yref:`u<L3Geom.u>` - :yref:`u0<L3Geom.u0>`); 1 means the true scale, 0 disables drawing local displacements; negative values are permissible."))
	);
	// clang-format on
};
REGISTER_SERIALIZABLE(Gl1_L3Geom);

struct Gl1_L6Geom : public Gl1_L3Geom {
	RENDERS(L6Geom);
	void go(const shared_ptr<IGeom>&, const shared_ptr<Interaction>&, const shared_ptr<Body>&, const shared_ptr<Body>&, bool) override;
	// clang-format off
	YADE_CLASS_BASE_DOC_STATICATTRS(Gl1_L6Geom,Gl1_L3Geom,"Render :yref:`L6Geom` geometry.",
		((Real,phiScale,1.,,"Scale local rotations (:yref:`phi<L6Geom.phi>` - :yref:`phi0<L6Geom.phi0>`). The default scale is to draw $\\pi$ rotation with length equal to minimum radius."))
	);
	// clang-format on
};
REGISTER_SERIALIZABLE(Gl1_L6Geom);
#endif

struct Ig2_Sphere_Sphere_L3Geom : public IGeomFunctor {
	virtual bool
	             go(const shared_ptr<Shape>&       s1,
	                const shared_ptr<Shape>&       s2,
	                const State&                   state1,
	                const State&                   state2,
	                const Vector3r&                shift2,
	                const bool&                    force,
	                const shared_ptr<Interaction>& I) override;
	virtual bool genericGo(
	        bool                           is6Dof,
	        const shared_ptr<Shape>&       s1,
	        const shared_ptr<Shape>&       s2,
	        const State&                   state1,
	        const State&                   state2,
	        const Vector3r&                shift2,
	        const bool&                    force,
	        const shared_ptr<Interaction>& I);
	// common code for {sphere,facet,wall}+sphere contacts
	// facet&wall will get separated if L3Geom subclass with exact branch vector is created
	void handleSpheresLikeContact(
	        const shared_ptr<Interaction>& I,
	        const State&                   state1,
	        const State&                   state2,
	        const Vector3r&                shift2,
	        bool                           is6Dof,
	        const Vector3r&                normal,
	        const Vector3r&                contPt,
	        Real                           uN,
	        Real                           r1,
	        Real                           r2);

	enum { APPROX_NO_MID_TRSF = 1, APPROX_NO_MID_NORMAL = 2, APPROX_NO_RENORM_MID_NORMAL = 4 };

	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS(Ig2_Sphere_Sphere_L3Geom,IGeomFunctor,"Incrementally compute :yref:`L3Geom` for contact of 2 spheres. Detailed documentation in py/\\_extraDocs.py",
		((bool,noRatch,true,,"See :yref:`Ig2_Sphere_Sphere_ScGeom.avoidGranularRatcheting`."))
		((Real,distFactor,1,,"Create interaction if spheres are not futher than *distFactor* \\*(r1+r2). If negative, zero normal deformation will be set to be the initial value (otherwise, the geometrical distance is the \'\'zero'' one)."))
		((int,trsfRenorm,100,,"How often to renormalize :yref:`trsf<L3Geom.trsf>`; if non-positive, never renormalized (simulation might be unstable)"))
		((int,approxMask,0,,"Selectively enable geometrical approximations (bitmask); add the values for approximations to be enabled.\n\n"
		"== ===============================================================\n"
		"1  use previous transformation to transform velocities (which are known at mid-steps), instead of mid-step transformation computed as quaternion slerp at t=0.5.\n"
		"2  do not take average (mid-step) normal when computing relative shear displacement, use previous value instead\n"
		"4  do not re-normalize average (mid-step) normal, if used.…\n"
		"== ===============================================================\n\n"
		"By default, the mask is zero, wherefore none of these approximations is used.\n"
		))
	);
	// clang-format on
	FUNCTOR2D(Sphere, Sphere);
	DEFINE_FUNCTOR_ORDER_2D(Sphere, Sphere);
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(Ig2_Sphere_Sphere_L3Geom);

struct Ig2_Wall_Sphere_L3Geom : public Ig2_Sphere_Sphere_L3Geom {
	virtual bool
	go(const shared_ptr<Shape>&       s1,
	   const shared_ptr<Shape>&       s2,
	   const State&                   state1,
	   const State&                   state2,
	   const Vector3r&                shift2,
	   const bool&                    force,
	   const shared_ptr<Interaction>& I) override;
	//virtual bool genericGo(bool is6Dof, const shared_ptr<Shape>& s1, const shared_ptr<Shape>& s2, const State& state1, const State& state2, const Vector3r& shift2, const bool& force, const shared_ptr<Interaction>& I);
	// clang-format off
	YADE_CLASS_BASE_DOC(Ig2_Wall_Sphere_L3Geom,Ig2_Sphere_Sphere_L3Geom,"Incrementally compute :yref:`L3Geom` for contact between :yref:`Wall` and :yref:`Sphere`. Uses attributes of :yref:`Ig2_Sphere_Sphere_L3Geom`.");
	// clang-format on
	FUNCTOR2D(Wall, Sphere);
	DEFINE_FUNCTOR_ORDER_2D(Wall, Sphere);
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(Ig2_Wall_Sphere_L3Geom);

#ifdef L3GEOM_SPHERESLIKE
struct Ig2_Facet_Sphere_L3Geom : public Ig2_Sphere_Sphere_L3Geom {
	// get point on segment A..B closest to P; algo: http://local.wasp.uwa.edu.au/~pbourke/geometry/pointline/
	static Vector3r getClosestSegmentPt(const Vector3r& P, const Vector3r& A, const Vector3r& B)
	{
		using math::max;
		using math::min;
		Vector3r BA = B - A;
		Real     u  = (P.dot(BA) - A.dot(BA)) / (BA.squaredNorm());
		return A + min((Real)1., max((Real)0., u)) * BA;
	}
	virtual bool
	go(const shared_ptr<Shape>&       s1,
	   const shared_ptr<Shape>&       s2,
	   const State&                   state1,
	   const State&                   state2,
	   const Vector3r&                shift2,
	   const bool&                    force,
	   const shared_ptr<Interaction>& I) override;
	// clang-format off
	YADE_CLASS_BASE_DOC(Ig2_Facet_Sphere_L3Geom,Ig2_Sphere_Sphere_L3Geom,"Incrementally compute :yref:`L3Geom` for contact between :yref:`Facet` and :yref:`Sphere`. Uses attributes of :yref:`Ig2_Sphere_Sphere_L3Geom`.");
	// clang-format on
	FUNCTOR2D(Facet, Sphere);
	DEFINE_FUNCTOR_ORDER_2D(Facet, Sphere);
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(Ig2_Facet_Sphere_L3Geom);
#endif

struct Ig2_Sphere_Sphere_L6Geom : public Ig2_Sphere_Sphere_L3Geom {
	virtual bool
	go(const shared_ptr<Shape>&       s1,
	   const shared_ptr<Shape>&       s2,
	   const State&                   state1,
	   const State&                   state2,
	   const Vector3r&                shift2,
	   const bool&                    force,
	   const shared_ptr<Interaction>& I) override;
	// clang-format off
	YADE_CLASS_BASE_DOC(Ig2_Sphere_Sphere_L6Geom,Ig2_Sphere_Sphere_L3Geom,"Incrementally compute :yref:`L6Geom` for contact of 2 spheres.");
	// clang-format on
	FUNCTOR2D(Sphere, Sphere);
	DEFINE_FUNCTOR_ORDER_2D(Sphere, Sphere);
};
REGISTER_SERIALIZABLE(Ig2_Sphere_Sphere_L6Geom);


struct Law2_L3Geom_FrictPhys_ElPerfPl : public LawFunctor {
	bool go(shared_ptr<IGeom>&, shared_ptr<IPhys>&, Interaction*) override;
	FUNCTOR2D(L3Geom, FrictPhys);
	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS(Law2_L3Geom_FrictPhys_ElPerfPl,LawFunctor,"Basic law for testing :yref:`L3Geom`; it bears no cohesion (unless *noBreak* is ``True``), and plastic slip obeys the Mohr-Coulomb criterion (unless *noSlip* is ``True``).",
		((bool,noBreak,false,,"Do not break contacts when particles separate."))
		((bool,noSlip,false,,"No plastic slipping."))
		((int,plastDissipIx,-1,(Attr::noSave|Attr::hidden),"Index of plastically dissipated energy"))
		((int,elastPotentialIx,-1,(Attr::hidden|Attr::noSave),"Index for elastic potential energy (with O.trackEnergy)"))
	);
	// clang-format on
};
REGISTER_SERIALIZABLE(Law2_L3Geom_FrictPhys_ElPerfPl);

struct Law2_L6Geom_FrictPhys_Linear : public Law2_L3Geom_FrictPhys_ElPerfPl {
	bool go(shared_ptr<IGeom>&, shared_ptr<IPhys>&, Interaction*) override;
	FUNCTOR2D(L6Geom, FrictPhys);
	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS(Law2_L6Geom_FrictPhys_Linear,Law2_L3Geom_FrictPhys_ElPerfPl,"Basic law for testing :yref:`L6Geom` -- linear in both normal and shear sense, without slip or breakage.",
		((Real,charLen,1,,"Characteristic length with the meaning of the stiffness ratios bending/shear and torsion/normal."))
	);
	// clang-format on
};
REGISTER_SERIALIZABLE(Law2_L6Geom_FrictPhys_Linear);

} // namespace yade
