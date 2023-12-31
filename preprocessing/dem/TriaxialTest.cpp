/*************************************************************************
*  Copyright (C) 2006 by Bruno Chareyre		                         *
*  bruno.chareyre@grenoble-inp.fr                                            *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/


#include <lib/high-precision/Constants.hpp>
#include <core/Aabb.hpp>
#include <core/Body.hpp>
#include <core/Dispatching.hpp>
#include <core/Interaction.hpp>
#include <core/InteractionLoop.hpp>
#include <core/Scene.hpp>
#include <pkg/common/Box.hpp>
#include <pkg/common/ElastMat.hpp>
#include <pkg/common/Facet.hpp>
#include <pkg/common/ForceResetter.hpp>
#include <pkg/common/GravityEngines.hpp>
#include <pkg/common/InsertionSortCollider.hpp>
#include <pkg/common/Sphere.hpp>
#include <pkg/common/Wall.hpp>
#include <pkg/dem/ElasticContactLaw.hpp>
#include <pkg/dem/FrictPhys.hpp>
#include <pkg/dem/GlobalStiffnessTimeStepper.hpp>
#include <pkg/dem/NewtonIntegrator.hpp>
#include <pkg/dem/TriaxialCompressionEngine.hpp>
#include <pkg/dem/TriaxialStateRecorder.hpp>
#include <pkg/dem/TriaxialStressController.hpp>
#include <preprocessing/dem/Shop.hpp>

#include <pkg/common/Bo1_Aabb.hpp>
#include <pkg/common/Wall.hpp>
#include <pkg/dem/Ig2_Box_Sphere_ScGeom.hpp>
#include <pkg/dem/Ig2_Facet_Sphere_ScGeom.hpp>
#include <pkg/dem/Ig2_Sphere_Sphere_ScGeom.hpp>

#include <boost/limits.hpp>
#include <boost/numeric/conversion/bounds.hpp>
#include <preprocessing/dem/SpherePack.hpp>
//#include<pkg/dem/MicroMacroAnalyser.hpp>

#include "TriaxialTest.hpp"

namespace yade { // Cannot have #include directive inside.

CREATE_LOGGER(TriaxialTest);
YADE_PLUGIN((TriaxialTest));

TriaxialTest::~TriaxialTest() { }

bool TriaxialTest::generate(string& message)
{
	message = "";
	if (facetWalls && wallWalls) { LOG_WARN("Turning TriaxialTest::facetWalls off, since wallWalls were selected as well."); }

	shared_ptr<Body> body;

	/* if _mean_radius is not given (i.e. <=0), then calculate it from box size;
	 * OTOH, if it is specified, scale the box preserving its ratio and lowerCorner so that the radius can be as requested
	 */
	Real       porosity = .8;
	SpherePack sphere_pack;
	if (importFilename == "") {
		Vector3r dimensions = upperCorner - lowerCorner;
		Real     volume     = dimensions.x() * dimensions.y() * dimensions.z();
		long     num;
		if (radiusMean <= 0) num = sphere_pack.makeCloud(lowerCorner, upperCorner, -1, radiusStdDev, numberOfGrains, false /*periodic?*/, porosity);
		else {
			bool fixedDims[3];
			fixedDims[0] = fixedBoxDims.find('x') != string::npos;
			fixedDims[1] = fixedBoxDims.find('y') != string::npos;
			fixedDims[2] = fixedBoxDims.find('z') != string::npos;
			int nScaled  = (3 - (int)fixedDims[0] + (int)fixedDims[1] + (int)fixedDims[2]);
			if (nScaled == 0)
				throw std::invalid_argument(
				        "At most 2 (not 3) axes can have fixed dimensions in fixedBoxDims if scaling for given radiusMean.");
			Real boxScaleFactor = radiusMean * pow((4 / 3.) * Mathr::PI * numberOfGrains / (volume * (1 - porosity)), 1. / nScaled);
			LOG_INFO("Mean radius value of " << radiusMean << " requested, scaling " << nScaled << " dimensions by " << boxScaleFactor);
			dimensions[0] *= fixedDims[0] ? 1. : boxScaleFactor;
			dimensions[1] *= fixedDims[1] ? 1. : boxScaleFactor;
			dimensions[2] *= fixedDims[2] ? 1. : boxScaleFactor;
			upperCorner = lowerCorner + dimensions;
			num         = sphere_pack.makeCloud(
                                lowerCorner, upperCorner, radiusMean, radiusStdDev, numberOfGrains, false, -1, vector<Real>(), vector<Real>(), false, seed);
		}
		message += "Generated a sample with " + boost::lexical_cast<string>(num) + " spheres inside box of dimensions: ("
		        + boost::lexical_cast<string>(upperCorner[0] - lowerCorner[0]) + "," + boost::lexical_cast<string>(upperCorner[1] - lowerCorner[1])
		        + "," + boost::lexical_cast<string>(upperCorner[2] - lowerCorner[2]) + ").";
	} else {
		if (radiusMean > 0) LOG_WARN("radiusMean ignored, since importFilename specified.");
		sphere_pack.fromFile(importFilename);
		sphere_pack.aabb(lowerCorner, upperCorner);
	}
	// setup scene here, since radiusMean is now at its true value (if it was negative)
	scene = shared_ptr<Scene>(new Scene);
	positionRootBody(scene);
	createActors(scene);

	if (thickness < 0) thickness = radiusMean;
	if (facetWalls || wallWalls) thickness = 0;
	if (!facetWalls && !wallWalls) {
		// bottom box
		Vector3r center   = Vector3r((lowerCorner[0] + upperCorner[0]) / 2, lowerCorner[1] - thickness / 2.0, (lowerCorner[2] + upperCorner[2]) / 2);
		Vector3r halfSize = Vector3r(
		        wallOversizeFactor * math::abs(lowerCorner[0] - upperCorner[0]) / 2 + thickness,
		        thickness / 2.0,
		        wallOversizeFactor * math::abs(lowerCorner[2] - upperCorner[2]) / 2 + thickness);
		createBox(body, center, halfSize, true);
		scene->bodies->insert(body);
		triaxialcompressionEngine->wall_bottom_id = body->getId();
		// top box
		center   = Vector3r((lowerCorner[0] + upperCorner[0]) / 2, upperCorner[1] + thickness / 2.0, (lowerCorner[2] + upperCorner[2]) / 2);
		halfSize = Vector3r(
		        wallOversizeFactor * math::abs(lowerCorner[0] - upperCorner[0]) / 2 + thickness,
		        thickness / 2.0,
		        wallOversizeFactor * math::abs(lowerCorner[2] - upperCorner[2]) / 2 + thickness);
		createBox(body, center, halfSize, true);
		scene->bodies->insert(body);
		triaxialcompressionEngine->wall_top_id = body->getId();
		// box 1
		center   = Vector3r(lowerCorner[0] - thickness / 2.0, (lowerCorner[1] + upperCorner[1]) / 2, (lowerCorner[2] + upperCorner[2]) / 2);
		halfSize = Vector3r(
		        thickness / 2.0,
		        wallOversizeFactor * math::abs(lowerCorner[1] - upperCorner[1]) / 2 + thickness,
		        wallOversizeFactor * math::abs(lowerCorner[2] - upperCorner[2]) / 2 + thickness);
		createBox(body, center, halfSize, true);
		scene->bodies->insert(body);
		triaxialcompressionEngine->wall_left_id = body->getId();
		// box 2
		center   = Vector3r(upperCorner[0] + thickness / 2.0, (lowerCorner[1] + upperCorner[1]) / 2, (lowerCorner[2] + upperCorner[2]) / 2);
		halfSize = Vector3r(
		        thickness / 2.0,
		        wallOversizeFactor * math::abs(lowerCorner[1] - upperCorner[1]) / 2 + thickness,
		        wallOversizeFactor * math::abs(lowerCorner[2] - upperCorner[2]) / 2 + thickness);

		createBox(body, center, halfSize, true);
		scene->bodies->insert(body);
		triaxialcompressionEngine->wall_right_id = body->getId();
		// box 3
		center   = Vector3r((lowerCorner[0] + upperCorner[0]) / 2, (lowerCorner[1] + upperCorner[1]) / 2, lowerCorner[2] - thickness / 2.0);
		halfSize = Vector3r(
		        wallOversizeFactor * math::abs(lowerCorner[0] - upperCorner[0]) / 2 + thickness,
		        wallOversizeFactor * math::abs(lowerCorner[1] - upperCorner[1]) / 2 + thickness,
		        thickness / 2.0);
		createBox(body, center, halfSize, true);
		scene->bodies->insert(body);
		triaxialcompressionEngine->wall_back_id = body->getId();
		// box 4
		center   = Vector3r((lowerCorner[0] + upperCorner[0]) / 2, (lowerCorner[1] + upperCorner[1]) / 2, upperCorner[2] + thickness / 2.0);
		halfSize = Vector3r(
		        wallOversizeFactor * math::abs(lowerCorner[0] - upperCorner[0]) / 2 + thickness,
		        wallOversizeFactor * math::abs(lowerCorner[1] - upperCorner[1]) / 2 + thickness,
		        thickness / 2.0);
		createBox(body, center, halfSize, true);
		scene->bodies->insert(body);
		triaxialcompressionEngine->wall_front_id = body->getId();
	}
	size_t imax = sphere_pack.pack.size();
	for (size_t i = 0; i < imax; i++) {
		const SpherePack::Sph& sp(sphere_pack.pack[i]);
		LOG_DEBUG("sphere (" << sp.c << " " << sp.r << ")");
		createSphere(body, sp.c, sp.r, false, true);
		scene->bodies->insert(body);
	}
	if (defaultDt < 0) {
		defaultDt                             = Shop::PWaveTimeStep(scene);
		scene->dt                             = defaultDt;
		globalStiffnessTimeStepper->defaultDt = defaultDt;
		LOG_INFO("Computed default (PWave) timestep " << defaultDt);
	}
	return true;
}

void TriaxialTest::createSphere(shared_ptr<Body>& body, Vector3r position, Real radius, bool /*big*/, bool /*dynamic*/)
{
	body            = shared_ptr<Body>(new Body);
	body->groupMask = 2;
	shared_ptr<Aabb>   aabb(new Aabb);
	shared_ptr<Sphere> iSphere(new Sphere);
	body->state->blockedDOFs = State::DOF_NONE;
	body->state->mass        = 4.0 / 3.0 * Mathr::PI * radius * radius * radius * density;
	body->state->inertia     = Vector3r(
                2.0 / 5.0 * body->state->mass * radius * radius,
                2.0 / 5.0 * body->state->mass * radius * radius,
                2.0 / 5.0 * body->state->mass * radius * radius);
	body->state->pos = position;
	shared_ptr<FrictMat> mat(new FrictMat);
	mat->young         = sphereYoungModulus;
	mat->poisson       = sphereKsDivKn;
	mat->frictionAngle = compactionFrictionDeg * Mathr::PI / 180.0;
	aabb->color        = Vector3r(0, 1, 0);
	iSphere->radius    = radius;
	//iSphere->color	= Vector3r(0.4,0.1,0.1);
	iSphere->color = Vector3r(math::unitRandom(), math::unitRandom(), math::unitRandom());
	iSphere->color.normalize();
	body->shape    = iSphere;
	body->bound    = aabb;
	body->material = mat;
}


void TriaxialTest::createBox(shared_ptr<Body>& body, Vector3r position, Vector3r extents, bool wire)
{
	body                     = shared_ptr<Body>(new Body);
	body->groupMask          = 2;
	body->state->blockedDOFs = State::DOF_ALL;
	shared_ptr<Aabb> aabb(new Aabb);
	aabb->color      = Vector3r(1, 0, 0);
	body->bound      = aabb;
	body->state->pos = position;
	shared_ptr<FrictMat> mat(new FrictMat);
	mat->young         = sphereYoungModulus;
	mat->poisson       = sphereKsDivKn;
	mat->frictionAngle = boxFrictionDeg * Mathr::PI / 180.0;
	body->material     = mat;
	if (!facetWalls && !wallWalls) {
		shared_ptr<Box> iBox(new Box);
		iBox->extents = extents;
		iBox->wire    = wire;
		iBox->color   = Vector3r(1, 1, 1);
		body->shape   = iBox;
	}
	// guess the orientation
	int ax0 = extents[0] == 0 ? 0 : (extents[1] == 0 ? 1 : 2);
	int ax1 = (ax0 + 1) % 3, ax2 = (ax0 + 2) % 3;
	if (facetWalls) {
		Vector3r corner = position - extents; // "lower right" corner, with 90 degrees
		Vector3r side1(Vector3r::Zero());
		side1[ax1] = 4 * extents[ax1];
		Vector3r side2(Vector3r::Zero());
		side2[ax2] = 4 * extents[ax2];
		Vector3r v[3];
		v[0]                  = corner;
		v[1]                  = corner + side1;
		v[2]                  = corner + side2;
		Vector3r          cog = Shop::inscribedCircleCenter(v[0], v[1], v[2]);
		shared_ptr<Facet> iFacet(new Facet);
		for (int i = 0; i < 3; i++) {
			iFacet->vertices[i] = v[i] - cog;
		}
		iFacet->color = Vector3r(1, 1, 1);
		body->shape   = iFacet;
	}
	if (wallWalls) {
		shared_ptr<Wall> wall(new Wall);
		wall->sense = 0; // interact from both sides, since unspecified here
		wall->axis  = ax0;
		body->shape = wall;
	}
}


void TriaxialTest::createActors(shared_ptr<Scene>& scene2)
{
	// declaration of ‘scene’ shadows a member of ‘yade::TriaxialTest’ [-Werror=shadow]
	shared_ptr<IGeomDispatcher> interactionGeometryDispatcher(new IGeomDispatcher);
	interactionGeometryDispatcher->add(new Ig2_Sphere_Sphere_ScGeom);
	interactionGeometryDispatcher->add(new Ig2_Facet_Sphere_ScGeom);
	interactionGeometryDispatcher->add(new Ig2_Box_Sphere_ScGeom);
	shared_ptr<IPhysDispatcher> interactionPhysicsDispatcher(new IPhysDispatcher);
	shared_ptr<IPhysFunctor>    ss(new Ip2_FrictMat_FrictMat_FrictPhys);
	interactionPhysicsDispatcher->add(ss);

	shared_ptr<GravityEngine> gravityCondition(new GravityEngine);
	gravityCondition->gravity = gravity;

	globalStiffnessTimeStepper                         = shared_ptr<GlobalStiffnessTimeStepper>(new GlobalStiffnessTimeStepper);
	globalStiffnessTimeStepper->timeStepUpdateInterval = timeStepUpdateInterval;
	globalStiffnessTimeStepper->defaultDt              = defaultDt;

	// moving walls to regulate the stress applied + compress when the packing is dense an stable
	//cerr << "triaxialcompressionEngine = shared_ptr<TriaxialCompressionEngine> (new TriaxialCompressionEngine);" << std::endl;
	triaxialcompressionEngine = shared_ptr<TriaxialCompressionEngine>(new TriaxialCompressionEngine);
	//This prevent the deprecation warning. In fact this preprocessor in itself is becoming deprecated
	triaxialcompressionEngine->warn                      = 1;
	triaxialcompressionEngine->stiffnessUpdateInterval   = wallStiffnessUpdateInterval; // = stiffness update interval
	triaxialcompressionEngine->radiusControlInterval     = radiusControlInterval;       // = stiffness update interval
	triaxialcompressionEngine->sigmaIsoCompaction        = sigmaIsoCompaction;
	triaxialcompressionEngine->sigmaLateralConfinement   = sigmaLateralConfinement;
	triaxialcompressionEngine->max_vel                   = maxWallVelocity;
	triaxialcompressionEngine->thickness                 = thickness;
	triaxialcompressionEngine->strainRate                = strainRate;
	triaxialcompressionEngine->StabilityCriterion        = StabilityCriterion;
	triaxialcompressionEngine->autoCompressionActivation = autoCompressionActivation;
	triaxialcompressionEngine->autoUnload                = autoUnload;
	triaxialcompressionEngine->autoStopSimulation        = autoStopSimulation;
	triaxialcompressionEngine->internalCompaction        = internalCompaction;
	triaxialcompressionEngine->maxMultiplier             = maxMultiplier;
	triaxialcompressionEngine->finalMaxMultiplier        = finalMaxMultiplier;
	triaxialcompressionEngine->Key                       = Key;
	triaxialcompressionEngine->noFiles                   = noFiles;
	triaxialcompressionEngine->frictionAngleDegree       = sphereFrictionDeg;
	triaxialcompressionEngine->fixedPoroCompaction       = false;
	triaxialcompressionEngine->fixedPorosity             = 1;
	// recording global stress
	if (recordIntervalIter > 0 && !noFiles) {
		triaxialStateRecorder             = shared_ptr<TriaxialStateRecorder>(new TriaxialStateRecorder);
		triaxialStateRecorder->file       = WallStressRecordFile + Key;
		triaxialStateRecorder->iterPeriod = recordIntervalIter;
	}
	scene2->engines.clear();
	scene2->engines.push_back(shared_ptr<Engine>(new ForceResetter));
	shared_ptr<InsertionSortCollider> collider(new InsertionSortCollider);
	scene2->engines.push_back(collider);
	collider->verletDist = .5 * radiusMean;
	collider->boundDispatcher->add(new Bo1_Sphere_Aabb);
	collider->boundDispatcher->add(new Bo1_Box_Aabb);
	collider->boundDispatcher->add(new Bo1_Facet_Aabb);
	collider->boundDispatcher->add(new Bo1_Wall_Aabb);

	shared_ptr<InteractionLoop> ids(new InteractionLoop);
	ids->geomDispatcher = interactionGeometryDispatcher;
	ids->physDispatcher = interactionPhysicsDispatcher;
	ids->lawDispatcher  = shared_ptr<LawDispatcher>(new LawDispatcher);
	shared_ptr<Law2_ScGeom_FrictPhys_CundallStrack> see(new Law2_ScGeom_FrictPhys_CundallStrack);
	ids->lawDispatcher->add(see);
	scene2->engines.push_back(ids);
	scene2->engines.push_back(globalStiffnessTimeStepper);
	scene2->engines.push_back(triaxialcompressionEngine);
	if (recordIntervalIter > 0 && !noFiles) scene2->engines.push_back(triaxialStateRecorder);

	shared_ptr<NewtonIntegrator> newton(new NewtonIntegrator);
	newton->damping = dampingForce;
	scene2->engines.push_back(newton);
}

void TriaxialTest::positionRootBody(shared_ptr<Scene>& /*scene*/) { }

} // namespace yade
