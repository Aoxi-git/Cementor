/*************************************************************************
*  Copyright (C) 2006 by luc Scholtes                                    *
*  luc.scholtes@hmg.inpg.fr                                              *
*  Copyright (C) 2008 by Bruno Chareyre                                  *
*  bruno.chareyre@grenoble-inp.fr                                            *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include "TriaxialStateRecorder.hpp"
#include <lib/high-precision/Constants.hpp>
// #include <pkg/dem/TriaxialCompressionEngine.hpp>
#include <core/Omega.hpp>
#include <core/Scene.hpp>
#include <pkg/common/Sphere.hpp>
#include <pkg/dem/FrictPhys.hpp>
#include <pkg/dem/ScGeom.hpp>
#include <pkg/dem/TriaxialStressController.hpp>
#include <preprocessing/dem/Shop.hpp>

namespace yade { // Cannot have #include directive inside.

CREATE_LOGGER(TriaxialStateRecorder);
TriaxialStateRecorder::~TriaxialStateRecorder() {};

void TriaxialStateRecorder::action()
{
	// at the beginning of the file; write column titles
	if (out.tellp() == 0) out << "iteration s11 s22 s33 e11 e22 e33 unb_force porosity kineticE" << endl;

	if (!triaxialStressController) {
		vector<shared_ptr<Engine>>::iterator itFirst = scene->engines.begin();
		vector<shared_ptr<Engine>>::iterator itLast  = scene->engines.end();
		for (; itFirst != itLast; ++itFirst) {
			if ((*itFirst)->getClassName() == "TriaxialCompressionEngine" || (*itFirst)->getClassName() == "ThreeDTriaxialEngine"
			    || (*itFirst)->getClassName() == "TriaxialStressController") {
				LOG_DEBUG("stress controller engine found");
				triaxialStressController = YADE_PTR_CAST<TriaxialStressController>(*itFirst);
				//triaxialCompressionEngine = shared_ptr<TriaxialCompressionEngine> (static_cast<TriaxialCompressionEngine*> ( (*itFirst).get()));
			}
		}
		if (!triaxialStressController) LOG_ERROR("stress controller engine NOT found");
	}
	if (!(scene->iter % triaxialStressController->computeStressStrainInterval == 0)) triaxialStressController->computeStressStrain();

	/// Compute porosity :
	Real Vs = 0;
	Real V  = (triaxialStressController->height) * (triaxialStressController->width) * (triaxialStressController->depth);

	for (const auto& b : *scene->bodies) {
		if (!(b) || b->isClump()) continue;
		if (b->isDynamic()) {
			//Sorry, the next string was commented, because it gave a Warning "unused variable v". Anton Gladky
			//const Vector3r& v = b->state->vel;
			Vs += 1.3333333 * Mathr::PI * pow(YADE_PTR_CAST<Sphere>(b->shape)->radius, 3);
		}
	}
	porosity = (V - Vs) / V;

	out << boost::lexical_cast<string>(scene->iter) << " "
	    << boost::lexical_cast<string>(triaxialStressController->stress[triaxialStressController->wall_right][0]) << " "
	    << boost::lexical_cast<string>(triaxialStressController->stress[triaxialStressController->wall_top][1]) << " "
	    << boost::lexical_cast<string>(triaxialStressController->stress[triaxialStressController->wall_front][2]) << " "
	    << boost::lexical_cast<string>(triaxialStressController->strain[0]) << " " << boost::lexical_cast<string>(triaxialStressController->strain[1])
	    << " " << boost::lexical_cast<string>(triaxialStressController->strain[2]) << " "
	    << boost::lexical_cast<string>(triaxialStressController->ComputeUnbalancedForce()) << " " << boost::lexical_cast<string>(porosity) << " "
	    << boost::lexical_cast<string>(Shop::kineticEnergy()) << endl;
}

YADE_PLUGIN((TriaxialStateRecorder));

} // namespace yade
