/*************************************************************************
*  Copyright (C) 2009 by Luc Sibille                                     *
*  luc.sibille@univ-nantes.fr                                            *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include "ThreeDTriaxialEngine.hpp"
#include <lib/high-precision/Constants.hpp>
#include <core/Omega.hpp>
#include <core/Scene.hpp>

#include <lib/base/Math.hpp>
#include <core/Interaction.hpp>
#include <pkg/common/ElastMat.hpp>
#include <pkg/common/Sphere.hpp>
#include <pkg/dem/FrictPhys.hpp>
#include <boost/lambda/lambda.hpp>
#include <preprocessing/dem/Shop.hpp>

namespace yade { // Cannot have #include directive inside.

class Ip2_CohFrictMat_CohFrictMat_CohFrictPhys;

CREATE_LOGGER(ThreeDTriaxialEngine);
YADE_PLUGIN((ThreeDTriaxialEngine));

ThreeDTriaxialEngine::~ThreeDTriaxialEngine() { }

void ThreeDTriaxialEngine::action()
{
	static int warn = 0;
	if (!warn++) LOG_WARN("This engine is deprecated, please switch to TriaxialStressController if you expect long term support.")
	if (firstRun) {
		LOG_INFO("First run, will initialize!");

		if (updateFrictionAngle) setContactProperties(frictionAngleDegree);

		height0 = height;
		depth0  = depth;
		width0  = width;

		if (stressControl_1) {
			wall_right_activated = true;
			wall_left_activated  = true; //are the right walls for direction 1?
		} else {
			wall_right_activated = false;
			wall_left_activated  = false;
		}

		if (stressControl_2) {
			wall_bottom_activated = true;
			wall_top_activated    = true;
		} else {
			wall_bottom_activated = false;
			wall_top_activated    = false;
		}

		if (stressControl_3) {
			wall_front_activated = true;
			wall_back_activated  = true; //are the right walls for direction 3?
		} else {
			wall_front_activated = false;
			wall_back_activated  = false;
		}

		//internalCompaction=false;  //is needed to avoid a control for internal compaction by the TriaxialStressController engine

		// 		isAxisymetric=false; //is needed to avoid a stress control according the parameter sigma_iso (but according to sigma1, sigma2 and sigma3)

		firstRun = false;
	}


	const Real& dt = scene->dt;

	if (!stressControl_1) // control in strain if wanted
	{
		if (currentStrainRate1 != strainRate1) currentStrainRate1 += (strainRate1 - currentStrainRate1) * (1 - strainDamping);

		State* p_left = Body::byId(wall_left_id, scene)->state.get();
		p_left->pos += 0.5 * currentStrainRate1 * width * translationAxisx * dt;
		State* p_right = Body::byId(wall_right_id, scene)->state.get();
		p_right->pos -= 0.5 * currentStrainRate1 * width * translationAxisx * dt;

	} else {
		if (currentStrainRate1 != strainRate1) currentStrainRate1 += (strainRate1 - currentStrainRate1) * (1 - strainDamping);
		max_vel1 = 0.5 * currentStrainRate1 * width;
	}


	if (!stressControl_2) // control in strain if wanted
	{
		if (currentStrainRate2 != strainRate2) currentStrainRate2 += (strainRate2 - currentStrainRate2) * (1 - strainDamping);

		State* p_bottom = Body::byId(wall_bottom_id, scene)->state.get();
		p_bottom->pos += 0.5 * currentStrainRate2 * height * translationAxisy * dt;
		State* p_top = Body::byId(wall_top_id, scene)->state.get();
		p_top->pos -= 0.5 * currentStrainRate2 * height * translationAxisy * dt;

	} else {
		if (currentStrainRate2 != strainRate2) currentStrainRate2 += (strainRate2 - currentStrainRate2) * (1 - strainDamping);
		max_vel2 = 0.5 * currentStrainRate2 * height;
	}


	if (!stressControl_3) // control in strain if wanted
	{
		if (currentStrainRate3 != strainRate3) currentStrainRate3 += (strainRate3 - currentStrainRate3) * (1 - strainDamping);


		State* p_back = Body::byId(wall_back_id, scene)->state.get();
		p_back->pos += 0.5 * currentStrainRate3 * depth * translationAxisz * dt;
		State* p_front = Body::byId(wall_front_id, scene)->state.get();
		p_front->pos -= 0.5 * currentStrainRate3 * depth * translationAxisz * dt;

	} else {
		if (currentStrainRate3 != strainRate3) currentStrainRate3 += (strainRate3 - currentStrainRate3) * (1 - strainDamping);
		max_vel3 = 0.5 * currentStrainRate3 * depth;
	}

	TriaxialStressController::action(); // this function is called to perform the external stress control or the internal compaction
}

void ThreeDTriaxialEngine::setContactProperties(Real frictionDegree)
{
	scene                             = Omega::instance().getScene().get();
	shared_ptr<BodyContainer>& bodies = scene->bodies;
	for (const auto& b : *scene->bodies) {
		if (b->isDynamic()) YADE_PTR_CAST<FrictMat>(b->material)->frictionAngle = frictionDegree * Mathr::PI / 180.0;
	}

	FOREACH(const shared_ptr<Interaction>& ii, *scene->interactions)
	{
		if (!ii->isReal()) continue;
		const shared_ptr<FrictMat>& sdec1 = YADE_PTR_CAST<FrictMat>((*bodies)[(Body::id_t)((ii)->getId1())]->material);
		const shared_ptr<FrictMat>& sdec2 = YADE_PTR_CAST<FrictMat>((*bodies)[(Body::id_t)((ii)->getId2())]->material);
		//FIXME - why dynamic_cast fails here?
		FrictPhys* contactPhysics              = YADE_CAST<FrictPhys*>((ii)->phys.get());
		Real       fa                          = sdec1->frictionAngle;
		Real       fb                          = sdec2->frictionAngle;
		contactPhysics->tangensOfFrictionAngle = math::tan(math::min(fa, fb));
	}
}

} // namespace yade
