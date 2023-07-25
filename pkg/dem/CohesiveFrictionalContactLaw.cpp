/*************************************************************************
*  Copyright (C) 2007 by Bruno Chareyre <bruno.chareyre@imag.fr>         *
*  Copyright (C) 2008 by Janek Kozicki <cosurgi@berlios.de>              *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include "CohesiveFrictionalContactLaw.hpp"
#include <core/Omega.hpp>
#include <core/Scene.hpp>
#include <pkg/dem/ScGeom.hpp>
#include <lib/base/LoggingUtils.hpp>

namespace yade { // Cannot have #include directive inside.

using math::max;
using math::min; // using inside .cpp file is ok.

YADE_PLUGIN((CohesiveFrictionalContactLaw)(Law2_ScGeom6D_CohFrictPhys_CohesionMoment)(CohFrictMat)(CohFrictPhys)(Ip2_CohFrictMat_CohFrictMat_CohFrictPhys));
CREATE_LOGGER(Law2_ScGeom6D_CohFrictPhys_CohesionMoment);

Real Law2_ScGeom6D_CohFrictPhys_CohesionMoment::getPlasticDissipation() const { return (Real)plasticDissipation; }
void Law2_ScGeom6D_CohFrictPhys_CohesionMoment::initPlasticDissipation(Real initVal)
{
	plasticDissipation.reset();
	plasticDissipation += initVal;
}

Real Law2_ScGeom6D_CohFrictPhys_CohesionMoment::normElastEnergy()
{
	Real normEnergy = 0;
	FOREACH(const shared_ptr<Interaction>& I, *scene->interactions)
	{
		if (!I->isReal()) continue;
		CohFrictPhys* phys = YADE_CAST<CohFrictPhys*>(I->phys.get());
		if (phys) { normEnergy += 0.5 * (phys->normalForce.squaredNorm() / phys->kn); }
	}
	return normEnergy;
}
Real Law2_ScGeom6D_CohFrictPhys_CohesionMoment::shearElastEnergy()
{
	Real shearEnergy = 0;
	FOREACH(const shared_ptr<Interaction>& I, *scene->interactions)
	{
		if (!I->isReal()) continue;
		CohFrictPhys* phys = YADE_CAST<CohFrictPhys*>(I->phys.get());
		if (phys) { shearEnergy += 0.5 * (phys->shearForce.squaredNorm() / phys->ks); }
	}
	return shearEnergy;
}

Real Law2_ScGeom6D_CohFrictPhys_CohesionMoment::bendingElastEnergy()
{
	Real bendingEnergy = 0;
	FOREACH(const shared_ptr<Interaction>& I, *scene->interactions)
	{
		if (!I->isReal()) continue;
		CohFrictPhys* phys = YADE_CAST<CohFrictPhys*>(I->phys.get());
		if (phys) { bendingEnergy += 0.5 * (phys->moment_bending.squaredNorm() / phys->kr); }
	}
	return bendingEnergy;
}

Real Law2_ScGeom6D_CohFrictPhys_CohesionMoment::twistElastEnergy()
{
	Real twistEnergy = 0;
	FOREACH(const shared_ptr<Interaction>& I, *scene->interactions)
	{
		if (!I->isReal()) continue;
		CohFrictPhys* phys = YADE_CAST<CohFrictPhys*>(I->phys.get());
		if (phys) { twistEnergy += 0.5 * (phys->moment_twist.squaredNorm() / phys->ktw); }
	}
	return twistEnergy;
}

Real Law2_ScGeom6D_CohFrictPhys_CohesionMoment::totalElastEnergy()
{
	Real totalEnergy = 0;
	FOREACH(const shared_ptr<Interaction>& I, *scene->interactions)
	{
		if (!I->isReal()) continue;
		CohFrictPhys* phys = YADE_CAST<CohFrictPhys*>(I->phys.get());
		if (phys) {
			totalEnergy += 0.5 * (phys->normalForce.squaredNorm() / phys->kn);
			totalEnergy += 0.5 * (phys->shearForce.squaredNorm() / phys->ks);
			totalEnergy += 0.5 * (phys->moment_bending.squaredNorm() / phys->kr);
			totalEnergy += 0.5 * (phys->moment_twist.squaredNorm() / phys->ktw);
		}
	}
	return totalEnergy;
}


void CohesiveFrictionalContactLaw::action()
{
	if (!functor) functor = shared_ptr<Law2_ScGeom6D_CohFrictPhys_CohesionMoment>(new Law2_ScGeom6D_CohFrictPhys_CohesionMoment);
	functor->always_use_moment_law = always_use_moment_law;
	functor->shear_creep           = shear_creep;
	functor->twist_creep           = twist_creep;
	functor->creep_viscosity       = creep_viscosity;
	functor->scene                 = scene;
	FOREACH(const shared_ptr<Interaction>& I, *scene->interactions)
	{
		if (!I->isReal()) continue;
		functor->go(I->geom, I->phys, I.get());
	}
}


bool Law2_ScGeom6D_CohFrictPhys_CohesionMoment::checkPlasticity(ScGeom6D* geom, CohFrictPhys* phys, Real& Fn, bool computeMoment)
{
	// if one force/torque breaks, skip computing the others since the contact properties have changed.
	// Update thresholds and compute again, with some recursive calls in section 1. Max recursion depth is 1.
	
	//  1. ____________________ Compare force components to their max value  _________________
	const bool brittle = phys->fragile and not phys->cohesionBroken;
	Real Fs = 0, maxFs = 0, scalarRoll=0, maxRoll=0, scalarTwist=0, maxTwist=0;
	
	// check tensile force magnitude and max plastic disp for normal direction
	bool breaks = (-Fn) > phys->normalAdhesion or (phys->unpMax >= 0 and -(geom->penetrationDepth + phys->normalAdhesion / phys->kn) > phys->unpMax);
	// no need to compute the rest if the contact breaks
	if ((brittle or phys->unpMax >= 0) and breaks) {
		phys->SetBreakingState(always_use_moment_law);
		if (geom->penetrationDepth < 0) return false;  // FIXME: for fragile behavior the dissipated energy is not correct (possibly), it should include the elastic energy lost in breakage 
		else return checkPlasticity(geom, phys, Fn, computeMoment);
	}
	// check max shear force
	if (phys->ks>0) 
	{
		Real maxFs = phys->shearAdhesion;
		if (!phys->cohesionDisablesFriction || phys->cohesionBroken) maxFs += Fn * phys->tangensOfFrictionAngle;
		maxFs = math::max((Real)0, maxFs);
		Fs = phys->shearForce.norm();
		breaks = breaks or Fs > maxFs;
		if (brittle and breaks) {
			phys->SetBreakingState(always_use_moment_law);
			if (geom->penetrationDepth < 0) return false;
			else return checkPlasticity(geom, phys, Fn, computeMoment);
		}
	}
	// check torques
	if (computeMoment)
	{
		// rolling torque
		if (phys->kr > 0 and (phys->maxRollPl >= 0. or phys->rollingAdhesion >=0)) {
			maxRoll = phys->rollingAdhesion;
			if (phys->cohesionBroken or !phys->cohesionDisablesFriction) maxRoll += phys->maxRollPl * Fn;
			maxRoll = math::max((Real)0, maxRoll);
			scalarRoll = phys->moment_bending.norm();
			breaks = breaks or (scalarRoll > maxRoll and phys->rollingAdhesion >0);
			if (brittle and breaks) {
				phys->SetBreakingState(always_use_moment_law);		
				if (geom->penetrationDepth < 0) return false;
				else return checkPlasticity(geom, phys, Fn, computeMoment);
			}
		}
		// twisting torque
		if (phys->ktw > 0 and (phys->maxTwistPl >= 0. or phys->twistingAdhesion >=0)) {
			maxTwist = phys->twistingAdhesion;
			if (phys->cohesionBroken or !phys->cohesionDisablesFriction) maxTwist += phys->maxTwistPl * Fn;
			maxTwist = math::max((Real)0, maxTwist);
			scalarTwist = phys->moment_twist.norm();
			breaks = breaks or (scalarTwist > maxTwist and phys->twistingAdhesion >0);
			if (brittle and breaks) {
				phys->SetBreakingState(always_use_moment_law);		
				if (geom->penetrationDepth < 0) return false;
				else return checkPlasticity(geom, phys, Fn, computeMoment);
			}
		}
	}	
	//  2. ____________________ Correct the forces and increment dissipated energy  _________________
	if ((-Fn) > phys->normalAdhesion) { //normal plasticity
		if (scene->trackEnergy || traceEnergy) {
			Real dissipated = ((1 / phys->kn) * (-Fn - phys->normalAdhesion)) *phys->normalAdhesion;
			plasticDissipation += dissipated;
			if (scene->trackEnergy) scene->energy->add(dissipated, "normalDissip", normalDissipIx, /*reset*/ false);
		}
		Fn  = -phys->normalAdhesion;
		phys->unp = geom->penetrationDepth + phys->normalAdhesion / phys->kn;
	}
	phys->normalForce = Fn * geom->normal;

	if (Fs > maxFs) { //Plasticity condition on shear force
		Real ratio = maxFs / Fs;
		Vector3r trialForce = phys->shearForce;
		phys->shearForce *= ratio;
		if (scene->trackEnergy || traceEnergy) {
			Real sheardissip = ((1 / phys->ks) * (trialForce - phys->shearForce)) /*plastic disp*/.dot(phys->shearForce) /*active force*/;
			if (sheardissip > 0) {
				plasticDissipation += sheardissip;
				if (scene->trackEnergy) scene->energy->add(sheardissip, "shearDissip", shearDissipIx, /*reset*/ false);
			}
		}
	}	
	if (computeMoment) {
		// limit rolling moment to the plastic value, if required
		if (scalarRoll > maxRoll) { // fix maximum rolling moment
			Real ratio = maxRoll / scalarRoll;
			phys->moment_bending *= ratio;
			if (scene->trackEnergy) {
				Real bendingdissip = ((1 / phys->kr) * (scalarRoll - maxRoll) * maxRoll) /*active force*/;
				if (bendingdissip > 0) scene->energy->add(bendingdissip, "bendingDissip", bendingDissipIx, /*reset*/ false);
			}
		}
		// limit twisting moment to the plastic value, if required
		if (scalarTwist > maxTwist) { // fix maximum rolling moment
			Real ratio = maxTwist / scalarTwist;
			phys->moment_twist *= ratio;
			if (scene->trackEnergy) {
				Real twistdissip = ((1 / phys->ktw) * (scalarTwist - maxTwist) * maxTwist) /*active force*/;
				if (twistdissip > 0) scene->energy->add(twistdissip, "twistDissip", twistDissipIx, /*reset*/ false);
			}
		}
	}
	return true;
}


void Law2_ScGeom6D_CohFrictPhys_CohesionMoment::setElasticForces(shared_ptr<IGeom>& ig, shared_ptr<IPhys>& ip, Interaction* I, bool computeMoment, Real& Fn, const Real& dt)
{
	State* de1        = Body::byId(I->getId1(), scene)->state.get();
	State* de2        = Body::byId(I->getId2(), scene)->state.get();
	ScGeom6D*     geom            = YADE_CAST<ScGeom6D*>(ig.get());
	CohFrictPhys* phys            = YADE_CAST<CohFrictPhys*>(ip.get());
	// NORMAL
	Fn = phys->kn * (geom->penetrationDepth - phys->unp);
	
	// SHEAR
// 	// update orientation
// 	Vector3r&       shearForce = geom->rotate(phys->shearForce);
	// increment
	const Vector3r& dus        = geom->shearIncrement();
	phys->shearForce -= phys->ks * dus;

	// TORQUES
	if (computeMoment) {
		if (!useIncrementalForm) {
			if (phys->ktw >0) phys->moment_twist   = (geom->getTwist() * phys->ktw) * geom->normal;
			else phys->moment_twist = Vector3r::Zero();
			if (phys->kr >0) phys->moment_bending = geom->getBending() * phys->kr;
			else phys->moment_bending = Vector3r::Zero();
		} else { // Use incremental formulation to compute moment_twis and moment_bending (no twist_creep is applied)
			Vector3r relAngVel = geom->getRelAngVel(de1, de2, dt);
			// *** Bending ***//
			if (phys->kr >0) {
				Vector3r relAngVelBend = relAngVel - geom->normal.dot(relAngVel) * geom->normal; // keep only the bending part
				// incremental formulation for the bending moment (as for the shear part)
				phys->moment_bending           = phys->moment_bending - phys->kr * relAngVelBend * dt;
			}
			else phys->moment_bending = Vector3r::Zero();
			// ----------------------------------------------------------------------------------------
			// *** Torsion ***//
			if (phys->ktw > 0) {
				Vector3r relAngVelTwist = geom->normal.dot(relAngVel) * geom->normal;
				Vector3r relRotTwist    = relAngVelTwist * dt; // component of relative rotation along n
				phys->moment_twist           = phys->moment_twist - phys->ktw * relRotTwist;
			}
			else phys->moment_twist = Vector3r::Zero();
		}
	}
}

bool Law2_ScGeom6D_CohFrictPhys_CohesionMoment::go(shared_ptr<IGeom>& ig, shared_ptr<IPhys>& ip, Interaction* contact)
{
	const Real&   dt              = scene->dt;
	const int&    id1             = contact->getId1();
	const int&    id2             = contact->getId2();
	State* de1        = Body::byId(id1, scene)->state.get();
	State* de2        = Body::byId(id2, scene)->state.get();
	ScGeom6D*     geom            = YADE_CAST<ScGeom6D*>(ig.get());
	CohFrictPhys* phys            = YADE_CAST<CohFrictPhys*>(ip.get());
	Vector3r&     shearForceFirst = phys->shearForce;
	bool computeMoment = phys->momentRotationLaw && (!phys->cohesionBroken || always_use_moment_law);

	if (consistencyCheck)
	{
		if (phys->maxRollPl >= 0. and phys->maxTwistPl < 0.) LOG_WARN_ONCE("rolling moment is bounded but twisting moment is not (interaction "<<id1<<"-"<<id2<<"), check etaRoll and etaTwist of the materials if it is not intentional.");
		if (phys->maxRollPl < 0. and phys->maxTwistPl >= 0.) LOG_WARN_ONCE("twisting moment is bounded but rolling moment is not (interaction "<<id1<<"-"<<id2<<"), check etaRoll and etaTwist of the materials if it is not intentional.");
		if	((phys->maxRollPl >= 0. or phys->rollingAdhesion >=0 or phys->maxTwistPl or phys->twistingAdhesion >= 0) and not useIncrementalForm) LOG_WARN_ONCE("If :yref:`Law2_ScGeom6D_CohFrictPhys_CohesionMoment::useIncrementalForm` is false, then plasticity will not be applied correctly to torques (the total formulation will not reproduce irreversibility).");
		if	((always_use_moment_law or phys->maxRollPl!=0 or phys->maxTwistPl != 0 or phys->ktw !=0 or phys->kr !=0) and not phys->momentRotationLaw) LOG_WARN_ONCE("Interaction "  <<id1<<"-"<<id2<< "has some parameters related to contact moment defined but it will be ignored because i.phys.momentRotationLaw is False. Make sure momentRotationLaw is True in the material class if moments are needed.");
		consistencyCheck = false;
	}
	
	if (contact->isFresh(scene)) shearForceFirst = Vector3r::Zero();
	
	//  ____________________ Creep current forces if necessary _________________
	if (shear_creep) shearForceFirst -= phys->ks * (shearForceFirst * dt / creep_viscosity);
	if (twist_creep ) {
		if (!useIncrementalForm) {
				Real        viscosity_twist     = creep_viscosity * math::pow((2 * math::min(geom->radius1, geom->radius2)), 2) / 16.0;
				Real        angle_twist_creeped = geom->getTwist() * (1 - dt / viscosity_twist);
				Quaternionr q_twist(AngleAxisr(geom->getTwist(), geom->normal));
				Quaternionr q_twist_creeped(AngleAxisr(angle_twist_creeped, geom->normal));
				Quaternionr q_twist_delta(q_twist_creeped * q_twist.conjugate());
				geom->twistCreep = geom->twistCreep * q_twist_delta;			
		} else LOG_WARN_ONCE("Law2_ScGeom6D_CohFrictPhys_CohesionMoment: no twis creep is included if the incremental formulation is used.");
	}
	
	//  ____________________  Update orientation ___________________________________
	geom->rotate(phys->shearForce);
	if (computeMoment and useIncrementalForm) {
		if (phys->kr >0) geom->rotate(phys->moment_bending);
		if (phys->ktw >0) geom->rotate(phys->moment_twist);
	}

	//  ____________________ Linear elasticity giving "trial" force/torques _________________
	
	Real Fn=0;
	setElasticForces(ig, ip, contact, computeMoment, Fn, dt);

	//  ____________________ Failure and plasticity _________________________________________

	bool alive = checkPlasticity(geom, phys, Fn, computeMoment);	
	
	// ____________________  Delete interaction if no contact  ______________________________
	if (not alive) return false;

	//  ____________________  Apply the force _______________________________________________
	applyForceAtContactPoint(
			-phys->normalForce - phys->shearForce,
			geom->contactPoint,
			id1,
			de1->se3.position,
			id2,
			de2->se3.position + (scene->isPeriodic ? scene->cell->intrShiftPos(contact->cellDist) : Vector3r::Zero()));
	
	if (computeMoment) {
		// Apply moments now
		Vector3r moment = phys->moment_twist + phys->moment_bending;
		scene->forces.addTorque(id1, -moment);
		scene->forces.addTorque(id2, moment);
	}
	return true;
}


void Ip2_CohFrictMat_CohFrictMat_CohFrictPhys::go(
        const shared_ptr<Material>& b1 // CohFrictMat
        ,
        const shared_ptr<Material>& b2 // CohFrictMat
        ,
        const shared_ptr<Interaction>& interaction)
{
	CohFrictMat* sdec1 = static_cast<CohFrictMat*>(b1.get());
	CohFrictMat* sdec2 = static_cast<CohFrictMat*>(b2.get());
	ScGeom6D*    geom  = YADE_CAST<ScGeom6D*>(interaction->geom.get());

	//Create cohesive interractions only once
	if (setCohesionNow && cohesionDefinitionIteration == -1) cohesionDefinitionIteration = scene->iter;
	if (setCohesionNow && cohesionDefinitionIteration != -1 && cohesionDefinitionIteration != scene->iter) {
		cohesionDefinitionIteration = -1;
		setCohesionNow              = 0;
	}

	if (geom) {
		const auto normalCohPreCalculated
		        = (normalCohesion) ? (*normalCohesion)(b1->id, b2->id) : math::min(sdec1->normalCohesion, sdec2->normalCohesion);
		const auto shearCohPreCalculated = (shearCohesion) ? (*shearCohesion)(b1->id, b2->id) : math::min(sdec1->shearCohesion, sdec2->shearCohesion);
		// the max stress in pure bending is 4*M/πr^3 = 4*M/(Ar) (for a circular cross-section), if it controls failure, max moment is (r/4)*normalAdhesion
		const auto rollingAdhPreCalculated  = (rollingCohesion) ? (*rollingCohesion)(b1->id, b2->id) : 0.25 * normalCohPreCalculated * pow(math::min(geom->radius2, geom->radius1), 3);
		// the max shear stress in pure twisting is 2*Mt/πr^3 = 2*Mt/(Ar) (for a circular cross-section), if it controls failure, max moment is (r/2)*shearAdhesion
		const auto twistingAdhPreCalculated  = (twistingCohesion) ? (*twistingCohesion)(b1->id, b2->id) : 0.5 * shearCohPreCalculated * pow(math::min(geom->radius2, geom->radius1), 3);

		if (!interaction->phys) {
			interaction->phys            = shared_ptr<CohFrictPhys>(new CohFrictPhys());
			CohFrictPhys* contactPhysics = YADE_CAST<CohFrictPhys*>(interaction->phys.get());
			Real          Ea             = sdec1->young;
			Real          Eb             = sdec2->young;
			Real          Va             = sdec1->poisson;
			Real          Vb             = sdec2->poisson;
			Real          Da             = geom->radius1;
			Real          Db             = geom->radius2;
			Real          fa             = sdec1->frictionAngle;
			Real          fb             = sdec2->frictionAngle;
			Real          Kn             = 2.0 * Ea * Da * Eb * Db / (Ea * Da + Eb * Db); //harmonic average of two stiffnesses
			Real          frictionAngle  = (!frictAngle) ? math::min(fa, fb) : (*frictAngle)(sdec1->id, sdec2->id, fa, fb);

			// harmonic average of alphas parameters
			Real AlphaKr, AlphaKtw;
			if (sdec1->alphaKr && sdec2->alphaKr) AlphaKr = 2.0 * sdec1->alphaKr * sdec2->alphaKr / (sdec1->alphaKr + sdec2->alphaKr);
			else
				AlphaKr = 0;
			if (sdec1->alphaKtw && sdec2->alphaKtw) AlphaKtw = 2.0 * sdec1->alphaKtw * sdec2->alphaKtw / (sdec1->alphaKtw + sdec2->alphaKtw);
			else
				AlphaKtw = 0;

			Real Ks;
			if (Va && Vb)
				Ks = 2.0 * Ea * Da * Va * Eb * Db * Vb
				        / (Ea * Da * Va + Eb * Db * Vb); //harmonic average of two stiffnesses with ks=V*kn for each sphere
			else
				Ks = 0;

			contactPhysics->kr                     = Da * Db * Ks * AlphaKr;
			contactPhysics->ktw                    = Da * Db * Ks * AlphaKtw;
			contactPhysics->tangensOfFrictionAngle = math::tan(frictionAngle);

			if ((setCohesionOnNewContacts || setCohesionNow) && sdec1->isCohesive && sdec2->isCohesive) {
				contactPhysics->cohesionBroken = false;
				contactPhysics->normalAdhesion = normalCohPreCalculated * pow(math::min(Db, Da), 2);
				contactPhysics->shearAdhesion  = shearCohPreCalculated * pow(math::min(Db, Da), 2);
				contactPhysics->rollingAdhesion  = rollingAdhPreCalculated * pow(math::min(geom->radius2, geom->radius1), 3);
				contactPhysics->twistingAdhesion  = twistingAdhPreCalculated * pow(math::min(geom->radius2, geom->radius1), 3);
				geom->initRotations(*(Body::byId(interaction->getId1(), scene)->state), *(Body::byId(interaction->getId2(), scene)->state));
				contactPhysics->fragile = (sdec1->fragile || sdec2->fragile);
			}
			contactPhysics->kn = Kn;
			contactPhysics->ks = Ks;

			contactPhysics->maxRollPl         = min(sdec1->etaRoll * Da, sdec2->etaRoll * Db);
			contactPhysics->maxTwistPl        = min(sdec1->etaTwist * Da, sdec2->etaTwist * Db);
			contactPhysics->momentRotationLaw = (sdec1->momentRotationLaw && sdec2->momentRotationLaw);
		} else { // !isNew, but if setCohesionNow, all contacts are initialized like if they were newly created
			CohFrictPhys* contactPhysics = YADE_CAST<CohFrictPhys*>(interaction->phys.get());
			if ((setCohesionNow && sdec1->isCohesive && sdec2->isCohesive) || contactPhysics->initCohesion) {
				contactPhysics->cohesionBroken = false;
				contactPhysics->normalAdhesion = normalCohPreCalculated * pow(math::min(geom->radius2, geom->radius1), 2);
				contactPhysics->shearAdhesion  = shearCohPreCalculated * pow(math::min(geom->radius2, geom->radius1), 2);
				contactPhysics->rollingAdhesion  = rollingAdhPreCalculated * pow(math::min(geom->radius2, geom->radius1), 3);
				contactPhysics->twistingAdhesion  = twistingAdhPreCalculated * pow(math::min(geom->radius2, geom->radius1), 3);
				geom->initRotations(*(Body::byId(interaction->getId1(), scene)->state), *(Body::byId(interaction->getId2(), scene)->state));
				contactPhysics->fragile      = (sdec1->fragile || sdec2->fragile);
				contactPhysics->initCohesion = false;
			}
		}
	}
};

} // namespace yade
