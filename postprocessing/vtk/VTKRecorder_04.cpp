#ifdef YADE_VTK

#include "VTKRecorder.hpp"
// https://codeyarns.com/2014/03/11/how-to-selectively-ignore-a-gcc-warning/
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpragmas"
#pragma GCC diagnostic ignored "-Wcomment"
#pragma GCC diagnostic ignored "-Wsuggest-override"
#include <lib/compatibility/VTKCompatibility.hpp> // fix InsertNextTupleValue → InsertNextTuple name change (and others in the future)

#include <vtkCellArray.h>
#include <vtkSmartPointer.h>

#ifdef YADE_MPI
#include <core/Subdomain.hpp>
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpragmas"
#pragma GCC diagnostic ignored "-Wcast-function-type"
#include <mpi.h>
#pragma GCC diagnostic pop
#endif

// Code that generates this warning, Note: we cannot do this trick in yade. If we have a warning in yade, we have to fix it! See also https://gitlab.com/yade-dev/trunk/merge_requests/73
// This method will work once g++ bug https://gcc.gnu.org/bugzilla/show_bug.cgi?id=53431#c34 is fixed.
#include <vtkHexahedron.h>
#include <vtkQuad.h>
#include <vtkTriangle.h>
#pragma GCC diagnostic pop

#include <core/Scene.hpp>
#include <pkg/common/Box.hpp>
#include <pkg/common/Facet.hpp>
#include <pkg/common/Sphere.hpp>
#include <pkg/dem/ConcretePM.hpp>
#include <pkg/dem/JointedCohesiveFrictionalPM.hpp>
#include <pkg/pfv/PartialSatClayEngine.hpp>
#ifdef YADE_LIQMIGRATION
#include <pkg/dem/ViscoelasticCapillarPM.hpp>
#endif
#include <boost/fusion/include/pair.hpp>
#include <boost/fusion/support/pair.hpp>
#include <boost/unordered_map.hpp>

namespace yade { // Cannot have #include directive inside.

#ifdef YADE_MASK_ARBITRARY
#define GET_MASK(b) b->groupMask.to_ulong()
#else
#define GET_MASK(b) b->groupMask
#endif

void VTKRecorder::action_04()
{
#ifdef YADE_MPI
	const auto& subD = YADE_PTR_CAST<Subdomain>(scene->subD);
	const auto& sz   = parallelMode ? subD->ids.size() : scene->bodies->size();
	for (unsigned bId = 0; bId != sz; ++bId) {
		const auto& b = parallelMode ? (*scene->bodies)[subD->ids[bId]] : (*scene->bodies)[bId];
#else
	for (const auto& b : *scene->bodies) {
#endif
		if (!b) continue;
		if (mask != 0 && !b->maskCompatible(mask)) continue;
		if (recActive[REC_SPHERES]) {
			const Sphere* sphere = dynamic_cast<Sphere*>(b->shape.get());
			if (sphere) {
				if (skipNondynamic && b->state->blockedDOFs == State::DOF_ALL) continue;
				vtkIdType pid[1];
				Vector3r  pos(scene->isPeriodic ? scene->cell->wrapShearedPt(b->state->pos) : b->state->pos);
				pid[0] = spheresPos->InsertNextPoint(pos);
				spheresCells->InsertNextCell(1, pid);
				radii->InsertNextValue(sphere->radius);
				if (recActive[REC_BSTRESS]) {
					const Matrix3r&                         bStress = bStresses[b->getId()];
					Eigen::SelfAdjointEigenSolver<Matrix3r> solver(
					        bStress); // bStress is probably not symmetric (= self-adjoint for real matrices), but the solver still works, considering only one half of bStress. Which is good since existence of (real) eigenvalues is not sure for not symmetric bStress..
					Matrix3r dirAll = solver.eigenvectors();
					Vector3r eigenVal
					        = solver.eigenvalues(); // cf http://eigen.tuxfamily.org/dox/classEigen_1_1SelfAdjointEigenSolver.html#a30caf3c3884a7f4a46b8ec94efd23c5e to be sure that eigenVal[i] * dirAll.col(i) = bStress * dirAll.col(i) and that eigenVal[0] <= eigenVal[1] <= eigenVal[2]
					spheresSigI->InsertNextValue(eigenVal[2]);
					spheresSigII->InsertNextValue(eigenVal[1]);
					spheresSigIII->InsertNextValue(eigenVal[0]);
					Vector3r dirI((Real)dirAll(0, 2), (Real)dirAll(1, 2), (Real)dirAll(2, 2));
					Vector3r dirII((Real)dirAll(0, 1), (Real)dirAll(1, 1), (Real)dirAll(2, 1));
					Vector3r dirIII((Real)dirAll(0, 0), (Real)dirAll(1, 0), (Real)dirAll(2, 0));
					spheresDirI->InsertNextTuple(dirI);
					spheresDirII->InsertNextTuple(dirII);
					spheresDirIII->InsertNextTuple(dirIII);
				}
				if (recActive[REC_ID]) spheresId->InsertNextValue(b->getId());
				if (recActive[REC_MASK]) spheresMask->InsertNextValue(GET_MASK(b));
				if (recActive[REC_MASS]) spheresMass->InsertNextValue(b->state->mass);
#ifdef YADE_MPI
				if (recActive[REC_SUBDOMAIN]) spheresSubdomain->InsertNextValue(b->subdomain);
#endif

#ifdef THERMAL
				if (recActive[REC_TEMP]) {
					auto* thState = b->state.get();
					spheresTemp->InsertNextValue(thState->temp);
				}
#endif
#ifdef PARTIALSAT
				if (recActive[REC_PARTIALSAT]) {
					PartialSatState* state = dynamic_cast<PartialSatState*>(b->state.get());
					spheresRadiiChange->InsertNextValue(state->radiiChange);
					spheresSuction->InsertNextValue(state->suction);
					spheresIncidentCells->InsertNextValue(state->lastIncidentCells);
				}
#endif
				if (recActive[REC_CLUMPID]) clumpId->InsertNextValue(b->clumpId);
				if (recActive[REC_COLORS]) {
					const Vector3r& color = sphere->color;
					spheresColors->InsertNextTuple(color);
				}
				if (recActive[REC_VELOCITY]) {
					Vector3r vel = Vector3r::Zero();
					if (scene->isPeriodic) { // Take care of cell deformation
						vel = scene->cell->bodyFluctuationVel(b->state->pos, b->state->vel, scene->cell->prevVelGrad)
						        + scene->cell->prevVelGrad * scene->cell->wrapShearedPt(b->state->pos);
					} else {
						vel = b->state->vel;
					}
					spheresLinVelVec->InsertNextTuple(vel);
					spheresLinVelLen->InsertNextValue(vel.norm());
					const Vector3r& angVel = b->state->angVel;
					spheresAngVelVec->InsertNextTuple(angVel);
					spheresAngVelLen->InsertNextValue(angVel.norm());
				}
				if (recActive[REC_STRESS]) {
					const Vector3r& stress = bodyStates[b->getId()].normStress;
					const Vector3r& shear  = bodyStates[b->getId()].shearStress;
					spheresNormalStressVec->InsertNextTuple(stress);
					spheresShearStressVec->InsertNextTuple(shear);
					spheresNormalStressNorm->InsertNextValue(stress.norm());
				}
				if (recActive[REC_LUBRICATION]) {
					const Matrix3r& ncs = NCStresses[b->getId()];
					const Matrix3r& scs = SCStresses[b->getId()];
					const Matrix3r& nls = NLStresses[b->getId()];
					const Matrix3r& sls = SLStresses[b->getId()];
					const Matrix3r& nps = NPStresses[b->getId()];

					spheresLubricationNormalContactStress->InsertNextTuple(ncs);
					spheresLubricationShearContactStress->InsertNextTuple(scs);
					spheresLubricationNormalLubricationStress->InsertNextTuple(nls);
					spheresLubricationShearLubricationStress->InsertNextTuple(sls);
					spheresLubricationNormalPotentialStress->InsertNextTuple(nps);
				}
				if (recActive[REC_FORCE]) {
					scene->forces.sync();
					const Vector3r& f  = scene->forces.getForce(b->getId());
					const Vector3r& t  = scene->forces.getTorque(b->getId());
					Real            fn = f.norm();
					Real            tn = t.norm();
					spheresForceLen->InsertNextValue(fn);
					spheresTorqueLen->InsertNextValue(tn);
					spheresForceVec->InsertNextTuple(f);
					spheresTorqueVec->InsertNextTuple(t);
				}

				if (recActive[REC_CPM]) {
					cpmDamage->InsertNextValue(YADE_PTR_CAST<CpmState>(b->state)->normDmg);
					const Matrix3r& ss = YADE_PTR_CAST<CpmState>(b->state)->stress;
					//Real s[3]={ss[0],ss[1],ss[2]};
					cpmStress->InsertNextTuple(ss);
				}

				if (recActive[REC_JCFPM]) {
					nbCracks->InsertNextValue(YADE_PTR_CAST<JCFpmState>(b->state)->nbBrokenBonds);
					jcfpmDamage->InsertNextValue(YADE_PTR_CAST<JCFpmState>(b->state)->damageIndex);
				}

				if (recActive[REC_COORDNUMBER]) { spheresCoordNumb->InsertNextValue(b->coordNumber()); }
#ifdef YADE_SPH
				if (recActive[REC_SPH]) {
					spheresRhoSPH->InsertNextValue(b->state->rho);
					spheresPressSPH->InsertNextValue(b->state->press);
					spheresCoordNumbSPH->InsertNextValue(b->coordNumber());
				}
#endif

#ifdef YADE_DEFORM
				if (recActive[REC_DEFORM]) {
					const Sphere* sphereDef = dynamic_cast<Sphere*>(b->shape.get());
					spheresRealRad->InsertNextValue(b->state->dR + sphereDef->radius);
				}
#endif

#ifdef YADE_LIQMIGRATION
				if (recActive[REC_LIQ]) {
					spheresLiqVol->InsertNextValue(b->state->Vf);
					const Real tmpVolIter = liqVolIterBody(b);
					spheresLiqVolIter->InsertNextValue(tmpVolIter);
					spheresLiqVolTotal->InsertNextValue(tmpVolIter + b->state->Vf);
				}
#endif
				if (recActive[REC_MATERIALID]) spheresMaterialId->InsertNextValue(b->material->id);
				continue;
			}
		} // end rec sphere.
		if (recActive[REC_FACETS]) {
			const Facet* facet = dynamic_cast<Facet*>(b->shape.get());
			if (facet) {
				Vector3r                     pos(scene->isPeriodic ? scene->cell->wrapShearedPt(b->state->pos) : b->state->pos);
				const vector<Vector3r>&      localPos   = facet->vertices;
				Matrix3r                     facetAxisT = b->state->ori.toRotationMatrix();
				vtkSmartPointer<vtkTriangle> tri        = vtkSmartPointer<vtkTriangle>::New();
				vtkIdType                    nbPoints   = facetsPos->GetNumberOfPoints();
				for (int i = 0; i < 3; ++i) {
					Vector3r globalPos = pos + facetAxisT * localPos[i];
					facetsPos->InsertNextPoint(globalPos);
					tri->GetPointIds()->SetId(i, nbPoints + i);
				}
				facetsCells->InsertNextCell(tri);
				if (recActive[REC_COLORS]) {
					const Vector3r& color = facet->color;
					facetsColors->InsertNextTuple(color);
				}
				if (recActive[REC_STRESS]) {
					const Vector3r& stress = bodyStates[b->getId()].normStress + bodyStates[b->getId()].shearStress;
					facetsStressVec->InsertNextTuple(stress);
					facetsStressLen->InsertNextValue(stress.norm());
				}
				if (recActive[REC_FORCE]) {
					scene->forces.sync();
					const Vector3r& f  = scene->forces.getForce(b->getId());
					const Vector3r& t  = scene->forces.getTorque(b->getId());
					Real            fn = f.norm();
					Real            tn = t.norm();
					facetsForceLen->InsertNextValue(fn);
					facetsTorqueLen->InsertNextValue(tn);
					facetsForceVec->InsertNextTuple(f);
					facetsTorqueVec->InsertNextTuple(t);
				}

				if (recActive[REC_MATERIALID]) facetsMaterialId->InsertNextValue(b->material->id);
				if (recActive[REC_MASK]) facetsMask->InsertNextValue(GET_MASK(b));
				if (recActive[REC_COORDNUMBER]) { facetsCoordNumb->InsertNextValue(b->coordNumber()); }
				continue;
			}
		}
		if (recActive[REC_BOXES]) {
			const Box* box = dynamic_cast<Box*>(b->shape.get());
			if (box) {
				Vector3r                 pos(scene->isPeriodic ? scene->cell->wrapShearedPt(b->state->pos) : b->state->pos);
				Quaternionr              ori(b->state->ori);
				Vector3r                 ext(box->extents);
				vtkSmartPointer<vtkQuad> boxes = vtkSmartPointer<vtkQuad>::New();
				Vector3r                 A     = Vector3r(-ext[0], -ext[1], -ext[2]);
				Vector3r                 B     = Vector3r(-ext[0], +ext[1], -ext[2]);
				Vector3r                 C     = Vector3r(+ext[0], +ext[1], -ext[2]);
				Vector3r                 D     = Vector3r(+ext[0], -ext[1], -ext[2]);

				Vector3r E = Vector3r(-ext[0], -ext[1], +ext[2]);
				Vector3r F = Vector3r(-ext[0], +ext[1], +ext[2]);
				Vector3r G = Vector3r(+ext[0], +ext[1], +ext[2]);
				Vector3r H = Vector3r(+ext[0], -ext[1], +ext[2]);

				A = pos + ori * A;
				B = pos + ori * B;
				C = pos + ori * C;
				D = pos + ori * D;
				E = pos + ori * E;
				F = pos + ori * F;
				G = pos + ori * G;
				H = pos + ori * H;

				addWallVTK(boxes, boxesPos, A, B, C, D);
				boxesCells->InsertNextCell(boxes);

				addWallVTK(boxes, boxesPos, E, H, G, F);
				boxesCells->InsertNextCell(boxes);

				addWallVTK(boxes, boxesPos, A, E, F, B);
				boxesCells->InsertNextCell(boxes);

				addWallVTK(boxes, boxesPos, G, H, D, C);
				boxesCells->InsertNextCell(boxes);

				addWallVTK(boxes, boxesPos, F, G, C, B);
				boxesCells->InsertNextCell(boxes);

				addWallVTK(boxes, boxesPos, D, H, E, A);
				boxesCells->InsertNextCell(boxes);

				for (int i = 0; i < 6; i++) {
					if (recActive[REC_COLORS]) {
						const Vector3r& color = box->color;
						boxesColors->InsertNextTuple(color);
					}
					if (recActive[REC_STRESS]) {
						const Vector3r& stress = bodyStates[b->getId()].normStress + bodyStates[b->getId()].shearStress;
						boxesStressVec->InsertNextTuple(stress);
						boxesStressLen->InsertNextValue(stress.norm());
					}
					if (recActive[REC_FORCE]) {
						scene->forces.sync();
						const Vector3r& f  = scene->forces.getForce(b->getId());
						const Vector3r& t  = scene->forces.getTorque(b->getId());
						Real            fn = f.norm();
						Real            tn = t.norm();
						boxesForceVec->InsertNextTuple(f);
						boxesTorqueVec->InsertNextTuple(t);
						boxesForceLen->InsertNextValue(fn);
						boxesTorqueLen->InsertNextValue(tn);
					}
					if (recActive[REC_MATERIALID]) boxesMaterialId->InsertNextValue(b->material->id);
					if (recActive[REC_MASK]) boxesMask->InsertNextValue(GET_MASK(b));
				}
				continue;
			}
		}
	} // end bodies loop.

#ifdef YADE_MPI
	if ((!parallelMode and recActive[REC_PERICELL]) or (scene->subdomain == 0 and recActive[REC_PERICELL]))
#else
	if (recActive[REC_PERICELL])
#endif
	{
		const Matrix3r& hSize = scene->cell->hSize;
		Vector3r        v0    = hSize * Vector3r(0, 0, 1);
		Vector3r        v1    = hSize * Vector3r(0, 1, 1);
		Vector3r        v2    = hSize * Vector3r(1, 1, 1);
		Vector3r        v3    = hSize * Vector3r(1, 0, 1);
		Vector3r        v4    = hSize * Vector3r(0, 0, 0);
		Vector3r        v5    = hSize * Vector3r(0, 1, 0);
		Vector3r        v6    = hSize * Vector3r(1, 1, 0);
		Vector3r        v7    = hSize * Vector3r(1, 0, 0);
		pericellPoints->InsertNextPoint(v0);
		pericellPoints->InsertNextPoint(v1);
		pericellPoints->InsertNextPoint(v2);
		pericellPoints->InsertNextPoint(v3);
		pericellPoints->InsertNextPoint(v4);
		pericellPoints->InsertNextPoint(v5);
		pericellPoints->InsertNextPoint(v6);
		pericellPoints->InsertNextPoint(v7);
		vtkSmartPointer<vtkHexahedron> h = vtkSmartPointer<vtkHexahedron>::New();
		vtkIdList*                     l = h->GetPointIds();
		for (int i = 0; i < 8; i++) {
			l->SetId(i, i);
		}
		pericellHexa->InsertNextCell(h);
	}
}
#undef GET_MASK

} // namespace yade

#endif /* YADE_VTK */
