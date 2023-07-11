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

bool Law2_ScGeom6D_CohFrictPhys_CohesionMoment::go(shared_ptr<IGeom>& ig, shared_ptr<IPhys>& ip, Interaction* contact)
{
	const Real&   dt              = scene->dt;
	const int&    id1             = contact->getId1();
	const int&    id2             = contact->getId2();
	ScGeom6D*     geom            = YADE_CAST<ScGeom6D*>(ig.get());
	CohFrictPhys* phys            = YADE_CAST<CohFrictPhys*>(ip.get());
	Vector3r&     shearForceFirst = phys->shearForce;

	if (consistencyCheck)
	{
		if (phys->maxRollPl >= 0. and phys->maxTwistPl < 0.) LOG_WARN_ONCE("rolling moment is bounded but twisting moment is not (interaction "<<id1<<"-"<<id2<<"), check etaRoll and etaTwist of the materials if it is not intentional.");
		if (phys->maxRollPl < 0. and phys->maxTwistPl >= 0.) LOG_WARN_ONCE("twisting moment is bounded but rolling moment is not (interaction "<<id1<<"-"<<id2<<"), check etaRoll and etaTwist of the materials if it is not intentional.");
		consistencyCheck = false;
	}
	
	if (contact->isFresh(scene)) shearForceFirst = Vector3r::Zero();
	const Real& un = geom->penetrationDepth;
	State* de1        = Body::byId(id1, scene)->state.get();
	State* de2        = Body::byId(id2, scene)->state.get();

	//  ____________________ 1. Linear elasticity giving "trial" force/torques _________________
	
	// NORMAL
	Real Fn = phys->kn * (un - phys->unp);
	
	// SHEAR
	// creep if necessecary
	if (shear_creep) shearForceFirst -= phys->ks * (shearForceFirst * dt / creep_viscosity);
	// update orientation
	Vector3r&       shearForce = geom->rotate(phys->shearForce);
	// increment
	const Vector3r& dus        = geom->shearIncrement();
	shearForce -= phys->ks * dus;

	// TORQUES
	bool computeMoment = phys->momentRotationLaw && (!phys->cohesionBroken || always_use_moment_law);
	if (computeMoment) {
		if (!useIncrementalForm) {
			if (twist_creep) {
				Real        viscosity_twist     = creep_viscosity * math::pow((2 * math::min(geom->radius1, geom->radius2)), 2) / 16.0;
				Real        angle_twist_creeped = geom->getTwist() * (1 - dt / viscosity_twist);
				Quaternionr q_twist(AngleAxisr(geom->getTwist(), geom->normal));
				Quaternionr q_twist_creeped(AngleAxisr(angle_twist_creeped, geom->normal));
				Quaternionr q_twist_delta(q_twist_creeped * q_twist.conjugate());
				geom->twistCreep = geom->twistCreep * q_twist_delta;
			}
			if (phys->ktw >0) phys->moment_twist   = (geom->getTwist() * phys->ktw) * geom->normal;
			else phys->moment_twist = Vector3r::Zero();
			if (phys->kr >0) phys->moment_bending = geom->getBending() * phys->kr;
			else phys->moment_bending = Vector3r::Zero();
		} else { // Use incremental formulation to compute moment_twis and moment_bending (no twist_creep is applied)
			if (twist_creep)
				LOG_WARN_ONCE("Law2_ScGeom6D_CohFrictPhys_CohesionMoment: no twis creep is included if the incremental formulation is used.");
			Vector3r relAngVel = geom->getRelAngVel(de1, de2, dt);
			// *** Bending ***//
			if (phys->kr >0) {
				Vector3r relAngVelBend = relAngVel - geom->normal.dot(relAngVel) * geom->normal; // keep only the bending part
				// incremental formulation for the bending moment (as for the shear part)
				phys->moment_bending           = geom->rotate(phys->moment_bending); // rotate moment vector (updated)
				phys->moment_bending           = phys->moment_bending - phys->kr * relAngVelBend * dt;
			}
			else phys->moment_bending = Vector3r::Zero();
			// ----------------------------------------------------------------------------------------
			// *** Torsion ***//
			if (phys->ktw > 0) {
				Vector3r relAngVelTwist = geom->normal.dot(relAngVel) * geom->normal;
				Vector3r relRotTwist    = relAngVelTwist * dt; // component of relative rotation along n
				// incremental formulation for the torsional moment
				phys->moment_twist           = geom->rotate(phys->moment_twist);             // rotate moment vector (updated)
				phys->moment_twist           = phys->moment_twist - phys->ktw * relRotTwist;
			}
			else phys->moment_twist = Vector3r::Zero();
		}
	}
	//  ____________________ 2. Failure and plasticity _________________
	
	

// 	void checkFailureCriteria(ScGeom6D* geom, CohFrictPhys* phys, Real Fn, bool computeMoment)
	
	// Check normal force first since other components depend on it
	if (phys->fragile && (-Fn) > phys->normalAdhesion) return false; // break due to tension, interaction will be reset
	
	if ((-Fn) > phys->normalAdhesion) { //normal plasticity
		Fn        = -phys->normalAdhesion;
		phys->unp = un + phys->normalAdhesion / phys->kn;
		if (phys->unpMax >= 0
			&& -phys->unp > phys->unpMax) // Actually unpMax should be defined as a function of the average particule sizes for instance
			return false;
	}
	phys->normalForce = Fn * geom->normal;
	
	Real Fs = phys->ks>0 ? phys->shearForce.norm() : 0;
	Real maxFs = 0;
	Real scalarRoll = 0;
	Real scalarTwist = 0;
	Real RollMax = 0;
	Real TwistMax = 0;
	if (computeMoment) {
		scalarRoll = phys->kr>0 ? phys->moment_bending.norm() : 0;
		scalarTwist = phys->ktw>0 ? phys->moment_twist.norm() : 0;
	}
	// if one component leads to failure, it can change the threshold for other components too, hence a "while" to recompute all when one fails
	bool allTested = false;
	while (not allTested) {
		// Max shear force		
		maxFs = phys->shearAdhesion;
		if (!phys->cohesionDisablesFriction || maxFs == 0) maxFs += Fn * phys->tangensOfFrictionAngle;
		maxFs = math::max((Real)0, maxFs);
		if (phys->fragile and not phys->cohesionBroken and Fs > maxFs) {
			phys->SetBreakingState(always_use_moment_law);
			continue;
		}
		if (computeMoment) {
			// Check rolling
			if (phys->kr > 0 and (phys->maxRollPl >= 0. or phys->rollingAdhesion >=0)) {
				RollMax = phys->rollingAdhesion;
				if (!phys->cohesionDisablesFriction or phys->cohesionBroken) RollMax += phys->maxRollPl * Fn;
				RollMax = math::max((Real)0, RollMax);
			}
			if (phys->fragile and not phys->cohesionBroken and scalarRoll > RollMax) {
				phys->SetBreakingState(always_use_moment_law);
				continue;
			}
			// Check twisting		
			if (phys->ktw > 0 and (phys->maxTwistPl >= 0. or phys->twistingAdhesion >=0)) {
				TwistMax = phys->twistingAdhesion;
				if (!phys->cohesionDisablesFriction or (phys->cohesionBroken and always_use_moment_law)) TwistMax += phys->maxTwistPl * Fn;
				TwistMax = math::max((Real)0, TwistMax);
			}
			if (phys->fragile and not phys->cohesionBroken and scalarTwist > TwistMax) {
				phys->SetBreakingState(always_use_moment_law);
				continue;
			}
		}
		allTested = true;
	}
	if (phys->cohesionBroken and Fn < 0) phys->normalForce = Vector3r::Zero(); // cancel the traction computed above if it's now broken
	
	if (Fs > maxFs) { //Plasticity condition on shear force
		maxFs               = maxFs / Fs;
		Vector3r trialForce = shearForce;
		shearForce *= maxFs;
		if (scene->trackEnergy || traceEnergy) {
			Real sheardissip = ((1 / phys->ks) * (trialForce - shearForce)) /*plastic disp*/.dot(shearForce) /*active force*/;
			if (sheardissip > 0) {
				plasticDissipation += sheardissip;
				if (scene->trackEnergy) scene->energy->add(sheardissip, "shearDissip", shearDissipIx, /*reset*/ false);
			}
		}
	}
	//Apply the force
	applyForceAtContactPoint(
			-phys->normalForce - shearForce,
			geom->contactPoint,
			id1,
			de1->se3.position,
			id2,
			de2->se3.position + (scene->isPeriodic ? scene->cell->intrShiftPos(contact->cellDist) : Vector3r::Zero()));

	if (computeMoment) {
		// limit rolling moment to the plastic value, if required
		if (!useIncrementalForm)
				LOG_WARN_ONCE("If :yref:`Law2_ScGeom6D_CohFrictPhys_CohesionMoment::useIncrementalForm` is false, then plasticity will not "
							"be applied correctly (the total formulation would not reproduce irreversibility).");
		if (scalarRoll > RollMax) { // fix maximum rolling moment
			Real ratio = RollMax / scalarRoll;
			phys->moment_bending *= ratio;
			if (scene->trackEnergy) {
				Real bendingdissip = ((1 / phys->kr) * (scalarRoll - RollMax) * RollMax) /*active force*/;
				if (bendingdissip > 0) scene->energy->add(bendingdissip, "bendingDissip", bendingDissipIx, /*reset*/ false);
			}
		}
		// limit twisting moment to the plastic value, if required
		if (!useIncrementalForm)
				LOG_WARN_ONCE("If :yref:`Law2_ScGeom6D_CohFrictPhys_CohesionMoment::useIncrementalForm` is false, then plasticity will not "
							"be applied correctly (the total formulation would not reproduce irreversibility).");
		if (scalarTwist > TwistMax) { // fix maximum rolling moment
			Real ratio = TwistMax / scalarTwist;
			phys->moment_twist *= ratio;
			if (scene->trackEnergy) {
				Real twistdissip = ((1 / phys->ktw) * (scalarTwist - TwistMax) * TwistMax) /*active force*/;
				if (twistdissip > 0) scene->energy->add(twistdissip, "twistDissip", twistDissipIx, /*reset*/ false);
			}
		}
		// Apply moments now
		Vector3r moment = phys->moment_twist + phys->moment_bending;
		scene->forces.addTorque(id1, -moment);
		scene->forces.addTorque(id2, moment);
	}
	/// Moment law END       ///
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
		const auto normalAdhPreCalculated
		        = (normalCohesion) ? (*normalCohesion)(b1->id, b2->id) : math::min(sdec1->normalCohesion, sdec2->normalCohesion);
		const auto shearAdhPreCalculated = (shearCohesion) ? (*shearCohesion)(b1->id, b2->id) : math::min(sdec1->shearCohesion, sdec2->shearCohesion);

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
				contactPhysics->normalAdhesion = normalAdhPreCalculated * pow(math::min(Db, Da), 2);
				contactPhysics->shearAdhesion  = shearAdhPreCalculated * pow(math::min(Db, Da), 2);
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
				contactPhysics->normalAdhesion = normalAdhPreCalculated * pow(math::min(geom->radius2, geom->radius1), 2);
				contactPhysics->shearAdhesion  = shearAdhPreCalculated * pow(math::min(geom->radius2, geom->radius1), 2);

				geom->initRotations(*(Body::byId(interaction->getId1(), scene)->state), *(Body::byId(interaction->getId2(), scene)->state));
				contactPhysics->fragile      = (sdec1->fragile || sdec2->fragile);
				contactPhysics->initCohesion = false;
			}
		}
	}
};

} // namespace yade
