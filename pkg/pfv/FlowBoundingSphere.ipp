/*************************************************************************
*  Copyright (C) 2009 by Emanuele Catalano <catalano@grenoble-inp.fr>    *
*  Copyright (C) 2009 by Bruno Chareyre <bruno.chareyre@grenoble-inp.fr> *
*  Copyright (C) 2012 by Donia Marzougui <donia.marzougui@grenoble-inp.fr>*
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/
#ifdef FLOW_ENGINE

// #define XVIEW
#include "FlowBoundingSphere.hpp" //include after #define XVIEW
#include "vector"
#include <lib/high-precision/Constants.hpp>
#include <assert.h>
#include <fstream>
#include <iostream>
#include <new>
#include <sys/stat.h>
#include <sys/types.h>
#include <utility>

#ifdef XVIEW
// #include "Vue3D.h" //FIXME implicit dependencies will look for this class (out of tree) even ifndef XVIEW
#endif

#ifdef YADE_OPENMP
#include <omp.h>
// #define GS_OPEN_MP //It should never be defined if Yade is not using openmp
#endif

// #define USE_FAST_MATH 1

namespace yade { // Cannot have #include directive inside.

namespace CGT {

	typedef vector<Real> VectorR;

	//! Use this factor, or minLength, to reduce max permeability values (see usage below))
	const Real minLength = 0.01; //percentage of mean rad

	//! Factors including the effect of 1/2 symmetry in hydraulic radii
	const Real multSym1 = 1 / pow(2, 0.25);
	const Real multSym2 = 1 / pow(4, 0.25);

#ifdef XVIEW
	Vue3D Vue1;
#endif
	template <class Tesselation> FlowBoundingSphere<Tesselation>::~FlowBoundingSphere() { }
	template <class Tesselation> FlowBoundingSphere<Tesselation>::FlowBoundingSphere()
	{
		xMin = 1000.0, xMax = -10000.0, yMin = 1000.0, yMax = -10000.0, zMin = 1000.0, zMax = -10000.0;
		currentTes = 0;
		nOfSpheres = 0;
		sectionArea = 0, Height = 0, vTotal = 0;
		vtkInfiniteVertices = 0, vtkInfiniteCells = 0;
		viscosity = 1;
		fluidBulkModulus = 0;
		tessBasedForce = true;
		for (int i = 0; i < 6; i++)
			boundsIds[i] = 0;
		minPermLength = 1e-6; // multiplier applied on throat radius to define a minimal throat length (escaping coincident points)
		slipBoundary = false; //no-slip/symmetry conditions on lateral boundaries
		tolerance = 1e-07;
		relax = 1.9;
		ks = 0;
		distanceCorrection = true;
		clampKValues = true;
		meanKStat = true;
		KOptFactor = 0;
		noCache = true;
		pressureChanged = false;
		computeAllCells
		        = true; //might be turned false IF the code is reorganized (we can make a separate function to compute unitForceVectors outside compute_Permeability) AND it really matters for CPU time
		debugOut = true;
		rAverage
		        = false; /** use the average between the effective radius (inscribed sphere in facet) and the equivalent (circle surface = facet fluid surface) **/
		OUTPUT_BOUDARIES_RADII = false;
		rAverage
		        = false; /** if true use the average between the effective radius (inscribed sphere in facet) and the equivalent (circle surface = facet fluid surface) **/
		// 	areaR2Permeability=true;
		permeabilityMap = false;
		computedOnce = false;
		minKdivKmean = 0.0001;
		maxKdivKmean = 100.;
		ompThreads = 1;
		errorCode = 0;
		pxpos = ppval = NULL;
	}

	template <class Tesselation> void FlowBoundingSphere<Tesselation>::resetNetwork()
	{
		T[currentTes].Clear();
		this->resetLinearSystem();
	}

	template <class Tesselation> void FlowBoundingSphere<Tesselation>::resetLinearSystem() { noCache = true; }

	template <class Tesselation> void FlowBoundingSphere<Tesselation>::averageRelativeCellVelocity()
	{
		if (noCache && T[!currentTes].Max_id() <= 0) return;
		RTriangulation&     Tri = lastSolution().Triangulation();
		Point               posAvFacet;
		int                 numCells = 0;
		Real                facetFlowRate = 0;
		FiniteCellsIterator cellEnd = Tri.finite_cells_end();
		for (FiniteCellsIterator cell = Tri.finite_cells_begin(); cell != cellEnd; cell++) {
			if (cell->info().isGhost or cell->info().isAlpha) continue;
			cell->info().averageVelocity() = CGAL::NULL_VECTOR;
			numCells++;
			Real totFlowRate = 0; //used to acount for influxes in elements where pressure is imposed
			for (int i = 0; i < 4; i++)
				if (!Tri.is_infinite(cell->neighbor(i))) {
					CVector Surfk = cell->info() - cell->neighbor(i)->info();
					Real    area = sqrt(Surfk.squared_length());
					Surfk = Surfk / area;
					CVector branch = cell->vertex(facetVertices[i][0])->point().point() - cell->info();
					posAvFacet = (Point)cell->info() + (branch * Surfk) * Surfk;
					facetFlowRate = (cell->info().kNorm())[i] * (cell->info().shiftedP() - cell->neighbor(i)->info().shiftedP());
					totFlowRate += facetFlowRate;
					cell->info().averageVelocity() = cell->info().averageVelocity() + (facetFlowRate) * (posAvFacet - CGAL::ORIGIN);
				}
			//This is the influx term
			if (cell->info().Pcondition)
				cell->info().averageVelocity() = cell->info().averageVelocity() - (totFlowRate) * ((Point)cell->info() - CGAL::ORIGIN);
			//now divide by volume
			if (cell->info().volume() == 0) cerr << "zero volume pore interrupting velocity calculation" << endl;
			else
				cell->info().averageVelocity() = cell->info().averageVelocity() / math::abs(cell->info().volume());
		}
	}


	template <class Tesselation> bool FlowBoundingSphere<Tesselation>::isOnSolid(Real X, Real Y, Real Z)
	{
		RTriangulation&     Tri = lastSolution().Triangulation();
		FiniteCellsIterator cellEnd = Tri.finiteCellsEnd();
		for (FiniteCellsIterator cell = Tri.finite_cells_begin(); cell != cellEnd; cell++) {
			for (int i = 0; i < 4; i++) {
				Real radius = sqrt(cell->vertex(i)->point().weight());
				if (X < (cell->vertex(i)->point().x() + radius) && X > (cell->vertex(i)->point().x() - radius)) {
					if (Y < (cell->vertex(i)->point().y() + radius) && Y > (cell->vertex(i)->point().y() - radius)) {
						if (Z < (cell->vertex(i)->point().z() + radius) && Z > (cell->vertex(i)->point().z() - radius)) { return true; }
					}
				}
			}
		}
		return false;
	}
	template <class Tesselation> void FlowBoundingSphere<Tesselation>::averageFluidVelocity()
	{
		averageRelativeCellVelocity();
		RTriangulation&        Tri = lastSolution().Triangulation();
		int                    numVertex = 0;
		FiniteVerticesIterator verticesEnd = Tri.finite_vertices_end();
		for (FiniteVerticesIterator vIt = Tri.finite_vertices_begin(); vIt != verticesEnd; vIt++) {
			numVertex++;
		}

		vector<Real>         volumes;
		vector<CGT::CVector> velocityVolumes;
		velocityVolumes.resize(numVertex);
		volumes.resize(numVertex);

		for (FiniteVerticesIterator vIt = Tri.finite_vertices_begin(); vIt != verticesEnd; vIt++) {
			velocityVolumes[vIt->info().id()] = CGAL::NULL_VECTOR;
			volumes[vIt->info().id()] = 0.f;
		}

		FiniteCellsIterator cellEnd = Tri.finiteCellsEnd();
		for (FiniteCellsIterator cell = Tri.finite_cells_begin(); cell != cellEnd; cell++) {
			if (cell->info().fictious() == 0) {
				for (int i = 0; i < 4; i++) {
					velocityVolumes[cell->vertex(i)->info().id()] = velocityVolumes[cell->vertex(i)->info().id()]
					        + cell->info().averageVelocity() * math::abs(cell->info().volume());
					volumes[cell->vertex(i)->info().id()] = volumes[cell->vertex(i)->info().id()] + math::abs(cell->info().volume());
				}
			}
		}

		std::ofstream fluid_vel("Velocity", std::ios::out);
		Real          Rx = (xMax - xMin) / 10;
		Real          Ry = (yMax - yMin) / 12;
		Real          Rz = (zMax - zMin) / 20;
		CellHandle    cellula;

		CVector velocity = CGAL::NULL_VECTOR;
		int     i = 0;
		for (Real X = xMin + Rx; X < xMax; X += Rx) {
			for (Real Y = yMin + Ry; Y < yMax; Y += Ry) {
				velocity = CGAL::NULL_VECTOR;
				i = 0;
				for (Real Z = zMin + Rz; Z < zMax; Z += Rz) {
					cellula = Tri.locate(Point(X, Y, Z));
					for (int y = 0; y < 4; y++) {
						if (!cellula->vertex(y)->info().isFictious) {
							velocity = velocity
							        + (velocityVolumes[cellula->vertex(y)->info().id()] / volumes[cellula->vertex(y)->info().id()]);
							i++;
						}
					}
				}
				velocity = velocity / i;
				fluid_vel << X << " " << Y << " " << velocity << endl;
			}
		}
	}
	template <class Tesselation> vector<Real> FlowBoundingSphere<Tesselation>::averageFluidVelocityOnSphere(unsigned int Id_sph)
	{ //FIXME: we are computing everything again for each other Id_sph...
		if (noCache && T[!currentTes].Max_id() <= 0) return vector<Real>(3, 0);
		averageRelativeCellVelocity();
		RTriangulation& Tri = lastSolution().Triangulation();
		Real            volumes;
		CGT::CVector    velocityVolumes;
		vector<Real>    result;
		result.resize(3);
		velocityVolumes = CGAL::NULL_VECTOR;
		volumes = 0.f;

		FiniteCellsIterator cellEnd = Tri.finite_cells_end();
		for (FiniteCellsIterator cell = Tri.finite_cells_begin(); cell != cellEnd; cell++) {
			if (cell->info().fictious() == 0) {
				for (unsigned int i = 0; i < 4; i++) {
					if (cell->vertex(i)->info().id() == Id_sph) {
						velocityVolumes = velocityVolumes + cell->info().averageVelocity() * math::abs(cell->info().volume());
						volumes = volumes + math::abs(cell->info().volume());
					}
				}
			}
		}

		for (int i = 0; i < 3; i++)
			result[i] += velocityVolumes[i] / volumes;
		return result;
	}

	template <class Tesselation> Real FlowBoundingSphere<Tesselation>::getPorePressure(Real X, Real Y, Real Z)
	{
		const RTriangulation& Tri = lastSolution().Triangulation();
		CellHandle            cell = Tri.locate(CGT::Sphere(X, Y, Z));
		return cell->info().p();
	}

	template <class Tesselation> Real FlowBoundingSphere<Tesselation>::getPoreTemperature(Real X, Real Y, Real Z)
	{
		const RTriangulation& Tri = lastSolution().Triangulation();
		CellHandle            cell = Tri.locate(CGT::Sphere(X, Y, Z));
		return cell->info().temp();
	}

	template <class Tesselation> int FlowBoundingSphere<Tesselation>::getCell(Real X, Real Y, Real Z)
	{
		RTriangulation& Tri = lastSolution().Triangulation();
		CellHandle      cell = Tri.locate(CGT::Sphere(X, Y, Z));
		return cell->info().id;
	}

	template <class Tesselation> void FlowBoundingSphere<Tesselation>::measurePressureProfile(Real WallUpy, Real WallDowny)
	{
		using math::max;
		using math::min;

		RTriangulation& Tri = lastSolution().Triangulation();
		CellHandle      permeameter;
		std::ofstream   capture("Pressure_profile", std::ios::app);
		int             intervals = 5;
		int             captures = 6;
		Real            Rz = (zMax - zMin) / intervals;
		Real            Ry = (WallUpy - WallDowny) / captures;
		Real            X = (xMax + xMin) / 2;
		Real            Y = WallDowny;
		Real            pressure = 0.f;
		int             cell = 0;
		for (int i = 0; i < captures; i++) {
			for (Real Z = min(zMin, zMax); Z <= max(zMin, zMax); Z += math::abs(Rz)) {
				permeameter = Tri.locate(CGT::Sphere(X, Y, Z));
				pressure += permeameter->info().p();
				cell++;
			}
			Y += Ry;
			capture << pressure / cell << endl;
		}
	}
	template <class Tesselation> Real FlowBoundingSphere<Tesselation>::averageSlicePressure(Real Y)
	{
		RTriangulation& Tri = lastSolution().Triangulation();
		Real            P_ave = 0.f;
		int             n = 0;
		Real            Ry = (yMax - yMin) / 30;
		Real            Rx = (xMax - xMin) / 30;
		Real            Rz = (zMax - zMin) / 30;
		for (Real X = xMin; X <= xMax + Ry / 10; X = X + Rx) {
			for (Real Z = zMin; Z <= zMax + Ry / 10; Z = Z + Rz) {
				P_ave += Tri.locate(CGT::Sphere(X, Y, Z))->info().p();
				n++;
			}
		}
		P_ave /= n;
		return P_ave;
	}
	template <class Tesselation> Real FlowBoundingSphere<Tesselation>::averagePressure()
	{
		RTriangulation& Tri = lastSolution().Triangulation();
		Real            P = 0.f, Ppond = 0.f, Vpond = 0.f;
		int             n = 0;
		for (FiniteCellsIterator cell = Tri.finite_cells_begin(); cell != Tri.finite_cells_end(); cell++) {
			P += cell->info().p();
			n++;
			Ppond += cell->info().p() * math::abs(cell->info().volume());
			Vpond += math::abs(cell->info().volume());
		}
		P /= n;
		Ppond /= Vpond;
		return Ppond;
	}


	template <class Tesselation> void FlowBoundingSphere<Tesselation>::computeFacetForcesWithCache(bool onlyCache)
	{
		RTriangulation& Tri = T[currentTes].Triangulation();
		CVector         nullVect(0, 0, 0);
		//reset forces
		if (!onlyCache)
			for (FiniteVerticesIterator v = Tri.finite_vertices_begin(); v != Tri.finite_vertices_end(); ++v)
				v->info().forces = nullVect;

#ifdef parallel_forces
		if (noCache) {
			perVertexUnitForce.clear();
			perVertexPressure.clear();
			perVertexUnitForce.resize(T[currentTes].maxId + 1);
			perVertexPressure.resize(T[currentTes].maxId + 1);
		}
#endif
		CellHandle   neighbourCell;
		VertexHandle mirrorVertex;
		CVector      tempVect;
		//FIXME : Ema, be carefull with this (noCache), it needs to be turned true after retriangulation
		if (noCache) {
			for (VCellIterator cellIt = T[currentTes].cellHandles.begin(); cellIt != T[currentTes].cellHandles.end(); cellIt++) {
				CellHandle& cell = *cellIt;
				//reset cache
				for (int k = 0; k < 4; k++)
					cell->info().unitForceVectors[k] = nullVect;

				for (int j = 0; j < 4; j++)
					if (!Tri.is_infinite(cell->neighbor(j))) {
						neighbourCell = cell->neighbor(j);
						const CVector& Surfk = cell->info().facetSurfaces[j];
						//FIXME : later compute that fluidSurf only once in hydraulicRadius, for now keep full surface not modified in cell->info for comparison with other forces schemes
						//The ratio void surface / facet surface
						//Area of the facet (i.e. the triangle)
						Real area = sqrt(Surfk.squared_length());
						if (area <= 0) cerr << "AREA <= 0!!" << endl;
						CVector                     facetNormal = Surfk / area;
						const std::vector<CVector>& crossSections = cell->info().facetSphereCrossSections;
						//This is the cross-sectional area of the throat
						CVector fluidSurfk = cell->info().facetSurfaces[j] * cell->info().facetFluidSurfacesRatio[j];
						/// handle fictious vertex since we can get the projected surface easily here
						if (cell->vertex(j)->info().isFictious) {
							//projection of facet on the boundary
							Real projSurf = math::abs(Surfk[boundary(cell->vertex(j)->info().id()).coordinate]);
							tempVect = -projSurf * boundary(cell->vertex(j)->info().id()).normal;
							cell->vertex(j)->info().forces = cell->vertex(j)->info().forces + tempVect * cell->info().p();
							//define the cached value for later use with cache*p
							cell->info().unitForceVectors[j] = cell->info().unitForceVectors[j] + tempVect;
						}
						/// Apply weighted forces f_k=sqRad_k/sumSqRad*f
						CVector facetUnitForce = -fluidSurfk * cell->info().solidSurfaces[j][3];
						CVector facetForce = cell->info().p() * facetUnitForce;
						for (int y = 0; y < 3; y++) {
							//1st the drag (viscous) force weighted by surface of spheres in the throat
							cell->vertex(facetVertices[j][y])->info().forces = cell->vertex(facetVertices[j][y])->info().forces
							        + facetForce * cell->info().solidSurfaces[j][y];
							//(add to cached value)
							cell->info().unitForceVectors[facetVertices[j][y]] = cell->info().unitForceVectors[facetVertices[j][y]]
							        + facetUnitForce * cell->info().solidSurfaces[j][y];
							//2nd the partial integral of pore pressure, which boils down to weighting by partial cross-sectional area
							//uncomment to get total force / comment to get only viscous forces (Bruno)
							if (!cell->vertex(facetVertices[j][y])->info().isFictious) {
								cell->vertex(facetVertices[j][y])->info().forces
								        = cell->vertex(facetVertices[j][y])->info().forces
								        - facetNormal * cell->info().p() * crossSections[j][y];
								//add to cached value
								cell->info().unitForceVectors[facetVertices[j][y]]
								        = cell->info().unitForceVectors[facetVertices[j][y]]
								        - facetNormal * crossSections[j][y];
							}
						}
#ifdef parallel_forces
						perVertexUnitForce[cell->vertex(j)->info().id()].push_back(&(cell->info().unitForceVectors[j]));
						perVertexPressure[cell->vertex(j)->info().id()].push_back(&(cell->info().p()));
#endif
					}
			}
			noCache = false; //cache should always be defined after execution of this function
			if (onlyCache) return;
		}

#ifndef parallel_forces
		else { //use cached values
			for (FiniteCellsIterator cell = Tri.finite_cells_begin(); cell != cellEnd; cell++)
				for (int yy = 0; yy < 4; yy++)
					cell->vertex(yy)->info().forces
					        = cell->vertex(yy)->info().forces + cell->info().unitForceVectors[yy] * cell->info().p();
		}

#else
#ifdef YADE_OPENMP
#pragma omp parallel for num_threads(ompThreads)
#endif
		for (int vn = 0; vn <= T[currentTes].maxId; vn++) {
			if (T[currentTes].vertexHandles[vn] == NULL) continue;
			VertexHandle& v = T[currentTes].vertexHandles[vn];
			const int&    id = v->info().id();
			CVector       tf(0, 0, 0);
			int           k = 0;
			for (vector<const Real*>::iterator c = perVertexPressure[id].begin(); c != perVertexPressure[id].end(); c++)
				tf = tf + (*(perVertexUnitForce[id][k++])) * (**c);
			v->info().forces = tf;
		}
#endif
		// 	}
		if (debugOut) {
			CVector totalForce = nullVect;
			for (FiniteVerticesIterator v = Tri.finite_vertices_begin(); v != Tri.finite_vertices_end(); ++v) {
				if (!v->info().isFictious) totalForce = totalForce + v->info().forces;
				else if (boundary(v->info().id()).flowCondition == 1)
					totalForce = totalForce + v->info().forces;
			}
			cout << "totalForce = " << totalForce << endl;
		}
	}

	template <class Tesselation>
	void FlowBoundingSphere<Tesselation>::applySinusoidalPressure(RTriangulation& Tri, Real amplitude, Real averagePressure, Real loadIntervals)
	{
		Real       step = 1 / loadIntervals;
		VectorCell tmpCells;
		tmpCells.resize(10000);
		VCellIterator cellsIt = tmpCells.begin();
		for (Real alpha = 0; alpha < 1.001; alpha += step) {
			VCellIterator cellsEnd = Tri.incident_cells(T[currentTes].vertexHandles[yMaxId], cellsIt);
			for (VCellIterator it = tmpCells.begin(); it != cellsEnd; it++) {
				if (!Tri.is_infinite(*it)) {
					Point&      p1 = (*it)->info();
					CellHandle& cell = *it;
					if (p1.x() < xMin) cell->info().p() = averagePressure + amplitude;
					else if (p1.x() > xMax)
						cell->info().p() = averagePressure - amplitude;
					else if (p1.x() > (xMin + alpha * (xMax - xMin)) && p1.x() < (xMin + (alpha + step) * (xMax - xMin)))
						cell->info().p() = averagePressure + (amplitude) * (cos(alpha * M_PI));
				}
			}
		}
	}

	template <class Tesselation> void FlowBoundingSphere<Tesselation>::applyUserDefinedPressure(RTriangulation& Tri, vector<Real>& xpos, vector<Real>& pval)
	{
		if (!(xpos.size() && xpos.size() == pval.size())) {
			cerr << "Wrong definition of boundary pressure, check input" << endl;
			return;
		}
		pxpos = &xpos;
		ppval = &pval;
		Real       dx = xpos[1] - xpos[0];
		Real       xinit = xpos[0];
		Real       xlast = xpos.back();
		VectorCell tmpCells;
		tmpCells.resize(10000);
		VCellIterator cellsEnd = Tri.incident_cells(T[currentTes].vertexHandles[yMaxId], tmpCells.begin());
		for (VCellIterator it = tmpCells.begin(); it != cellsEnd; it++) {
			if (Tri.is_infinite(*it)) continue;
			Point&      p1 = (*it)->info();
			CellHandle& cell = *it;
			if (p1.x() < xinit || p1.x() > xlast) cerr << "udef pressure: cell out of range" << endl;
			else {
				Real frac, intg;
				frac = modf((p1.x() - xinit) / dx, &intg);
				cell->info().p() = pval[size_t(intg)] * (1 - frac) + pval[size_t(intg) + 1] * frac;
			}
		}
	}

	template <class Tesselation> CVector FlowBoundingSphere<Tesselation>::cellBarycenter(CellHandle& cell)
	{
		CVector center(0, 0, 0);
		for (int k = 0; k < 4; k++)
			center = center + 0.25 * (cell->vertex(k)->point().point() - CGAL::ORIGIN);
		return center;
	}

	template <class Tesselation> void FlowBoundingSphere<Tesselation>::interpolate(Tesselation& Tes, Tesselation& NewTes)
	{
		CellHandle      oldCell;
		RTriangulation& Tri = Tes.Triangulation();
#ifdef YADE_OPENMP
		const long size = NewTes.cellHandles.size();
#pragma omp parallel for num_threads(ompThreads > 0 ? ompThreads : 1)
		for (long i = 0; i < size; i++) {
			CellHandle& newCell = NewTes.cellHandles[i];
#else
		for (typename VectorCell::iterator cellIt = NewTes.cellHandles.begin(); cellIt != NewTes.cellHandles.end(); cellIt++) {
			CellHandle& newCell = *cellIt;
#endif
			if (newCell->info().isGhost) continue;
			CVector center(0, 0, 0);
			if (newCell->info().fictious() == 0)
				for (int k = 0; k < 4; k++)
					center = center + 0.25 * (Tes.vertex(newCell->vertex(k)->info().id())->point().point() - CGAL::ORIGIN);
			else {
				Real boundPos = 0;
				int  coord = 0;
				for (int k = 0; k < 4; k++) {
					if (!newCell->vertex(k)->info().isFictious)
						center = center
						        + (1. / (4. - newCell->info().fictious()))
						                * (Tes.vertex(newCell->vertex(k)->info().id())->point().point() - CGAL::ORIGIN);
				}
				for (int k = 0; k < 4; k++) {
					if (newCell->vertex(k)->info().isFictious) {
						coord = boundary(newCell->vertex(k)->info().id()).coordinate;
						boundPos = boundary(newCell->vertex(k)->info().id()).p[coord];
						center = CVector(
						        coord == 0 ? boundPos : center[0],
						        coord == 1 ? boundPos : center[1],
						        coord == 2 ? boundPos : center[2]);
					}
				}
			}
			oldCell = Tri.locate(CGT::Sphere(center[0], center[1], center[2]));
			if (!newCell->info().Pcondition) newCell->info().getInfo(oldCell->info());
			if (!newCell->info().Tcondition && thermalEngine) newCell->info().temp() = oldCell->info().temp();
			newCell->info().blocked = oldCell->info().blocked;
			// if (oldCell->info().isCavity) newCell->info().p()=oldCell->info().p(); // needed?
		}
	}

	template <class Tesselation> Real FlowBoundingSphere<Tesselation>::checkSphereFacetOverlap(const Sphere& v0, const Sphere& v1, const Sphere& v2)
	{
		//First, check that v0 projection fall between v1 and v2...
		Real dist = (v0.point() - v1.point()) * (v2.point() - v1.point());
		if (dist < 0) return 0;
		Real v1v2 = (v2.point() - v1.point()).squared_length();
		if (dist > v1v2) return 0;
		//... then, check distance
		Real m = (cross_product(v0.point() - v1.point(), v2.point() - v1.point())).squared_length() / v1v2;
		if (m < v0.weight()) {
			Real d = 2 * sqrt((v0.weight() - m));
			Real teta = 2 * acos(sqrt(m / v0.weight()));
			return 0.5 * (teta * v0.weight() - d * sqrt(m)); //this is S0, we use crossSection to avoid computing an "asin"
			                                                 // 		return crossSection-m*d;
		} else
			return 0;
	}

	template <class Tesselation> void FlowBoundingSphere<Tesselation>::setBlocked(CellHandle& cell)
	{
		RTriangulation& Tri = T[currentTes].Triangulation();
		if (cell->info().Pcondition) cell->info().p() = 0;
		else
			blockedCells.push_back(cell);
		for (int j = 0; j < 4; j++) {
			(cell->info().kNorm())[j] = 0;
			(cell->neighbor(j)->info().kNorm())[Tri.mirror_index(cell, j)] = 0;
		}
	}


	template <class Tesselation> void FlowBoundingSphere<Tesselation>::computePermeability()
	{
		using math::max;
		using math::min;

		if (debugOut) cout << "----Computing_Permeability------" << endl;
		RTriangulation& Tri = T[currentTes].Triangulation();
		VSolidTot = 0, Vtotalissimo = 0, vPoral = 0, sSolidTot = 0, vTotalPorosity = 0, vPoralPorosity = 0;
		FiniteCellsIterator cellEnd = Tri.finite_cells_end();

		CellHandle neighbourCell;

		Real k = 0, distance = 0, radius = 0;
		int  surfneg = 0;
		int  NEG = 0, POS = 0, pass = 0;

		bool ref = Tri.finite_cells_begin()->info().isvisited;
		Real meanK = 0, STDEV = 0, meanRadius = 0, meanDistance = 0;
		Real infiniteK = 1e10;

		for (VCellIterator cellIt = T[currentTes].cellHandles.begin(); cellIt != T[currentTes].cellHandles.end(); cellIt++) {
			CellHandle& cell = *cellIt;
			if (cell->info().blocked) {
				setBlocked(cell);
				cell->info().isvisited = !ref;
			}
			Point& p1 = cell->info();
			for (int j = 0; j < 4; j++) {
				neighbourCell = cell->neighbor(j);
				Point& p2 = neighbourCell->info();
				if (!Tri.is_infinite(neighbourCell) && (neighbourCell->info().isvisited == ref || computeAllCells)) {
					//compute and store the area of sphere-facet intersections for later use
					VertexHandle W[3];
					for (int kk = 0; kk < 3; kk++) {
						W[kk] = cell->vertex(facetVertices[j][kk]);
					}
					Sphere& v0 = W[0]->point();
					Sphere& v1 = W[1]->point();
					Sphere& v2 = W[2]->point();
					cell->info().facetSphereCrossSections[j] = CVector(
					        W[0]->info().isFictious ? 0
					                                : 0.5 * v0.weight()
					                        * acos((v1.point() - v0.point()) * (v2.point() - v0.point())
					                               / sqrt((v1.point() - v0.point()).squared_length()
					                                      * (v2.point() - v0.point()).squared_length())),
					        W[1]->info().isFictious ? 0
					                                : 0.5 * v1.weight()
					                        * acos((v0.point() - v1.point()) * (v2.point() - v1.point())
					                               / sqrt((v1.point() - v0.point()).squared_length()
					                                      * (v2.point() - v1.point()).squared_length())),
					        W[2]->info().isFictious ? 0
					                                : 0.5 * v2.weight()
					                        * acos((v0.point() - v2.point()) * (v1.point() - v2.point())
					                               / sqrt((v1.point() - v2.point()).squared_length()
					                                      * (v2.point() - v0.point()).squared_length())));
					//FIXME: it should be possible to skip completely blocked cells, currently the problem is it segfault for undefined areas
					// 				if (cell->info().blocked) continue;//We don't need permeability for blocked cells, it will be set to zero anyway
					pass += 1;
					CVector l = p1 - p2;
					distance = sqrt(l.squared_length());
					if (!rAverage) radius = 2 * computeHydraulicRadius(cell, j);
					else
						radius = (computeEffectiveRadius(cell, j) + computeEquivalentRadius(cell, j)) * 0.5;
					if (radius < 0) NEG++;
					else
						POS++;
					if (radius == 0) { cout << "INS-INS PROBLEM!!!!!!!" << endl; }
					Real fluidArea = 0;
					if (distance != 0) {
						if (minPermLength > 0 && distanceCorrection) distance = max(minPermLength * radius, distance);
						const CVector& Surfk = cell->info().facetSurfaces[j];
						Real           area = sqrt(Surfk.squared_length());
						const CVector& crossSections = cell->info().facetSphereCrossSections[j];
						Real           S0 = 0;
						S0 = checkSphereFacetOverlap(v0, v1, v2);
						if (S0 == 0) S0 = checkSphereFacetOverlap(v1, v2, v0);
						if (S0 == 0) S0 = checkSphereFacetOverlap(v2, v0, v1);
						//take absolute value, since in rare cases the surface can be negative (overlaping spheres)
						fluidArea = math::abs(area - crossSections[0] - crossSections[1] - crossSections[2] + S0);
						cell->info().facetFluidSurfacesRatio[j] = fluidArea / area;
						// kFactor<0 means we replace Poiseuille by Darcy localy, yielding a particle size-independent bulk conductivity
						if (kFactor > 0) cell->info().kNorm()[j] = kFactor * (fluidArea * pow(radius, 2)) / (8 * viscosity * distance);
						else
							cell->info().kNorm()[j] = -kFactor * area / distance;
						if (tempDependentViscosity && kFactor < 0
						    && thermalEngine) { //negative kFactor takes on direct permeability value instead of k/viscosity as above
							const Real avgTemp = (cell->info().temp() + neighbourCell->info().temp()) / 2.;
							const Real viscosityWithTemp
							        = -1.11962e-5 * avgTemp + 1.24084e-3; // linear approx of viscosity between 20-70 degC
							cell->info().kNorm()[j] = -kFactor * area / (viscosityWithTemp * distance);
						}
						if (cell->info().isCavity && neighbourCell->info().isCavity) {
							cell->info().kNorm()[j] = cavityFactor; // arbitrarily high conductivity for cavity neighbors
						}
						meanDistance += distance;
						meanRadius += radius;
						if (!cell->info().isCavity) meanK += (cell->info().kNorm())[j];

						if (!neighbourCell->info().isGhost)
							(neighbourCell->info().kNorm())[Tri.mirror_index(cell, j)] = (cell->info().kNorm())[j];
						if (k < 0 && debugOut) {
							surfneg += 1;
							cout << "__ k<0 __" << k << " "
							     << " fluidArea " << fluidArea << " area " << area << " " << crossSections[0] << " "
							     << crossSections[1] << " " << crossSections[2] << " " << W[0]->info().id() << " "
							     << W[1]->info().id() << " " << W[2]->info().id() << " " << p1 << " " << p2 << " test " << endl;
						}
					} else {
						cout << "infinite K1!" << endl;
						k = infiniteK;
						cell->info().kNorm()[j] = infiniteK;
					} //Will be corrected in the next loop
					if (!neighbourCell->info().isGhost)
						(neighbourCell->info().kNorm())[Tri.mirror_index(cell, j)] = (cell->info().kNorm())[j];
				}
			}
			cell->info().isvisited = !ref;
		}
		if (debugOut) cout << "surfneg est " << surfneg << endl;
		meanK /= pass;
		meanRadius /= pass;
		meanDistance /= pass;
		Real globalK;
		if (kFactor > 0)
			globalK = kFactor * meanDistance * vPoral
			        / (sSolidTot * 8. * viscosity); //An approximate value of macroscopic permeability, for clamping local values below
		else
			globalK = meanK;
		if (debugOut) {
			cout << "PassCompK = " << pass << endl;
			cout << "meanK = " << meanK << endl;
			cout << "globalK = " << globalK << endl;
			cout << "maxKdivKmean*globalK = " << maxKdivKmean * globalK << endl;
			cout << "minKdivKmean*globalK = " << minKdivKmean * globalK << endl;
			cout << "meanTubesRadius = " << meanRadius << endl;
			cout << "meanDistance = " << meanDistance << endl;
		}
		ref = Tri.finite_cells_begin()->info().isvisited;
		pass = 0;

		if (clampKValues)
			for (VCellIterator cellIt = T[currentTes].cellHandles.begin(); cellIt != T[currentTes].cellHandles.end(); cellIt++) {
				CellHandle& cell = *cellIt;
				for (int j = 0; j < 4; j++) {
					neighbourCell = cell->neighbor(j);
					if (!Tri.is_infinite(neighbourCell) && neighbourCell->info().isvisited == ref) {
						pass++;
						(cell->info().kNorm())[j] = max(minKdivKmean * globalK, min((cell->info().kNorm())[j], maxKdivKmean * globalK));
						(neighbourCell->info().kNorm())[Tri.mirror_index(cell, j)] = (cell->info().kNorm())[j];
					}
				}
			}
		if (debugOut) cout << "PassKcorrect = " << pass << endl;
		if (debugOut) cout << "POS = " << POS << " NEG = " << NEG << " pass = " << pass << endl;

		// A loop to compute the standard deviation of the local K distribution, and use it to include/exclude K values higher then (meanK +/- K_opt_factor*STDEV)
		if (meanKStat) {
			std::ofstream k_opt_file("k_stdev.txt", std::ios::out);
			ref = Tri.finite_cells_begin()->info().isvisited;
			pass = 0;
			for (FiniteCellsIterator cell = Tri.finite_cells_begin(); cell != cellEnd; cell++) {
				for (int j = 0; j < 4; j++) {
					neighbourCell = cell->neighbor(j);
					if (!Tri.is_infinite(neighbourCell) && neighbourCell->info().isvisited == ref) {
						pass++;
						STDEV += pow(((cell->info().kNorm())[j] - meanK), 2);
					}
				}
				cell->info().isvisited = !ref;
			}
			STDEV = sqrt(STDEV / pass);
			if (debugOut) cout << "PassSTDEV = " << pass << endl << "STATISTIC K" << endl;
			Real k_min = 0, k_max = meanK + KOptFactor * STDEV;
			cout << "Kmoy = " << meanK << " Standard Deviation = " << STDEV << endl << "kmin = " << k_min << " kmax = " << k_max << endl;
			ref = Tri.finite_cells_begin()->info().isvisited;
			pass = 0;
			for (FiniteCellsIterator cell = Tri.finite_cells_begin(); cell != cellEnd; cell++) {
				for (int j = 0; j < 4; j++) {
					neighbourCell = cell->neighbor(j);
					if (!Tri.is_infinite(neighbourCell) && neighbourCell->info().isvisited == ref) {
						pass += 1;
						if ((cell->info().kNorm())[j] > k_max) {
							(cell->info().kNorm())[j] = k_max;
							(neighbourCell->info().kNorm())[Tri.mirror_index(cell, j)] = (cell->info().kNorm())[j];
						}
						k_opt_file << KOptFactor << " " << (cell->info().kNorm())[j] << endl;
					}
				}
				cell->info().isvisited = !ref;
			}
			if (debugOut) cout << "PassKopt = " << pass << endl;
		}
		if (debugOut) {
			FiniteVerticesIterator verticesEnd = Tri.finite_vertices_end();
			Real                   Vgrains = 0;
			int                    grains = 0;
			for (FiniteVerticesIterator vIt = Tri.finite_vertices_begin(); vIt != verticesEnd; vIt++) {
				if (!vIt->info().isFictious && !vIt->info().isGhost) {
					grains += 1;
					Vgrains += 1.33333333 * M_PI * pow(vIt->point().weight(), 1.5);
				}
			}
			cout << grains << "grains - "
			     << "vTotal = " << vTotal << " Vgrains = " << Vgrains << " vPoral1 = " << (vTotal - Vgrains) << endl;
			cout << "Vtotalissimo = " << Vtotalissimo / 2 << " VSolidTot = " << VSolidTot / 2 << " vPoral2 = " << vPoral / 2
			     << " sSolidTot = " << sSolidTot << endl
			     << endl;
			if (!rAverage) cout << "------Hydraulic Radius is used for permeability computation------" << endl << endl;
			else
				cout << "------Average Radius is used for permeability computation------" << endl << endl;
			cout << "-----computed_Permeability-----" << endl;
		}
	}

	template <class Tesselation> vector<Real> FlowBoundingSphere<Tesselation>::getConstrictions()
	{
		RTriangulation&            Tri = T[currentTes].Triangulation();
		vector<Real>               constrictions;
		CellHandle                 neighbourCell;
		const FiniteCellsIterator& cellEnd = Tri.finite_cells_end();
		for (FiniteCellsIterator cell = Tri.finite_cells_begin(); cell != cellEnd; cell++) {
			if (cell->info().isGhost) continue; // retain only the cells with barycenter in the (0,0,0) period
			for (int j = 0; j < 4; j++) {
				neighbourCell = cell->neighbor(j);
				if (cell->info().id < neighbourCell->info().id) constrictions.push_back(computeEffectiveRadius(cell, j));
			}
		}
		return constrictions;
	}

	template <class Tesselation> vector<Constriction> FlowBoundingSphere<Tesselation>::getConstrictionsFull()
	{
		RTriangulation&            Tri = T[currentTes].Triangulation();
		vector<Constriction>       constrictions;
		CellHandle                 neighbourCell;
		const FiniteCellsIterator& cellEnd = Tri.finite_cells_end();
		for (FiniteCellsIterator cell = Tri.finite_cells_begin(); cell != cellEnd; cell++) {
			if (cell->info().isGhost) continue; // retain only the cells with barycenter in the (0,0,0) period
			for (int j = 0; j < 4; j++) {
				neighbourCell = cell->neighbor(j);
				if (cell->info().id < neighbourCell->info().id) {
					vector<Real>   rn;
					const CVector& normal = cell->info().facetSurfaces[j];
					if (!normal[0] && !normal[1] && !normal[2]) continue;
					rn.push_back(computeEffectiveRadius(cell, j));
					rn.push_back(normal[0]);
					rn.push_back(normal[1]);
					rn.push_back(normal[2]);
					Constriction cons(pair<int, int>(cell->info().id, neighbourCell->info().id), rn);
					constrictions.push_back(cons);
				}
			}
		}
		return constrictions;
	}

	template <class Tesselation> Real FlowBoundingSphere<Tesselation>::computeEffectiveRadius(CellHandle cell, int j)
	{
		RTriangulation& Tri = T[currentTes].Triangulation();
		if (Tri.is_infinite(cell->neighbor(j))) return 0;

		Point pos[3]; //spheres pos
		Real  r[3];   //spheres radius
		for (int i = 0; i < 3; i++) {
			pos[i] = cell->vertex(facetVertices[j][i])->point().point();
			r[i] = sqrt(cell->vertex(facetVertices[j][i])->point().weight());
		}

		Real reff = computeEffectiveRadiusByPosRadius(pos[0], r[0], pos[1], r[1], pos[2], r[2]);
		if (reff < 0) return 0; //happens very rarely, with bounding spheres most probably
		//if the facet involves one ore more bounding sphere, we return R with a minus sign
		if (cell->vertex(facetVertices[j][2])->info().isFictious || cell->vertex(facetVertices[j][1])->info().isFictious
		    || cell->vertex(facetVertices[j][2])->info().isFictious)
			return -reff;
		else
			return reff;
	}
	////compute inscribed radius independently by position and radius
	template <class Tesselation>
	Real FlowBoundingSphere<Tesselation>::computeEffectiveRadiusByPosRadius(
	        const Point& posA, const Real& rA, const Point& posB, const Real& rB, const Point& posC, const Real& rC)
	{
		CVector B = posB - posA;
		CVector x = B / sqrt(B.squared_length());
		CVector C = posC - posA;
		CVector z = CGAL::cross_product(x, C);
		CVector y = CGAL::cross_product(x, z);
		y = y / sqrt(y.squared_length());

		Real b1[2];
		b1[0] = B * x;
		b1[1] = B * y;
		Real c1[2];
		c1[0] = C * x;
		c1[1] = C * y;

		Real A = ((pow(rA, 2)) * (1 - c1[0] / b1[0]) + ((pow(rB, 2) * c1[0]) / b1[0]) - pow(rC, 2) + pow(c1[0], 2) + pow(c1[1], 2)
		          - ((pow(b1[0], 2) + pow(b1[1], 2)) * c1[0] / b1[0]))
		        / (2 * c1[1] - 2 * b1[1] * c1[0] / b1[0]);
		Real BB = (rA - rC - ((rA - rB) * c1[0] / b1[0])) / (c1[1] - b1[1] * c1[0] / b1[0]);
		Real CC = (pow(rA, 2) - pow(rB, 2) + pow(b1[0], 2) + pow(b1[1], 2)) / (2 * b1[0]);
		Real D = (rA - rB) / b1[0];
		Real E = b1[1] / b1[0];
		Real F = pow(CC, 2) + pow(E, 2) * pow(A, 2) - 2 * CC * E * A;

		Real c = -F - pow(A, 2) + pow(rA, 2);
		Real b = 2 * rA - 2 * (D - BB * E) * (CC - E * A) - 2 * A * BB;
		Real a = 1 - pow((D - BB * E), 2) - pow(BB, 2);

		if ((pow(b, 2) - 4 * a * c) < 0) { cout << "NEGATIVE DETERMINANT" << endl; }
		Real reff = (-b + sqrt(pow(b, 2) - 4 * a * c)) / (2 * a);
		return reff;
	}

	template <class Tesselation> Real FlowBoundingSphere<Tesselation>::computeEquivalentRadius(CellHandle cell, int j)
	{
		Real fluidSurf = sqrt(cell->info().facetSurfaces[j].squared_length()) * cell->info().facetFluidSurfacesRatio[j];
		return sqrt(fluidSurf / M_PI);
	}
	template <class Tesselation> Real FlowBoundingSphere<Tesselation>::computeHydraulicRadius(CellHandle cell, int j)
	{
		RTriangulation& Tri = T[currentTes].Triangulation();
		if (Tri.is_infinite(cell->neighbor(j))) return 0;
		Real Vpore = this->volumePoreVoronoiFraction(cell, j);
		Real Ssolid = this->surfaceSolidThroat(cell, j, slipBoundary, /*reuse the same facet data*/ true);

		//handle symmetry (tested ok)
		if (slipBoundary && facetNFictious > 0) {
			//! Include a multiplier so that permeability will be K/2 or K/4 in symmetry conditions
			Real mult = facetNFictious == 1 ? multSym1 : multSym2;
			return Vpore / Ssolid * mult;
		}
		return Vpore / Ssolid;
	}
	template <class Tesselation> void FlowBoundingSphere<Tesselation>::initializePressure(Real pZero)
	{
		RTriangulation&     Tri = T[currentTes].Triangulation();
		FiniteCellsIterator cellEnd = Tri.finite_cells_end();

		for (FiniteCellsIterator cell = Tri.finite_cells_begin(); cell != cellEnd; cell++) {
			if (!cell->info().Pcondition) cell->info().p() = pZero;
			cell->info().dv() = 0;
		}
		// cuboid bcs
		if (alphaBound < 0) {
			for (int bound = 0; bound < 6; bound++) {
				int& id = *boundsIds[bound];
				boundingCells[bound].clear();
				if (id < 0) continue;
				Boundary& bi = boundary(id);
				if (!bi.flowCondition) {
					VectorCell tmpCells;
					tmpCells.resize(10000);
					VCellIterator cells_it = tmpCells.begin();
					VCellIterator cells_end = Tri.incident_cells(T[currentTes].vertexHandles[id], cells_it);
					for (VCellIterator it = tmpCells.begin(); it != cells_end; it++) {
						(*it)->info().p() = bi.value;
						(*it)->info().Pcondition = true;
						boundingCells[bound].push_back(*it);
					}
				}
				boundingCells[bound].shrink_to_fit();
			}
		} else {
			// identify cells incident to alpha vertices and set BCs
			Tesselation& Tes = T[currentTes];
			const long   sizeCells = Tes.cellHandles.size();
			//#ifdef YADE_OPENMP
			//#pragma omp parallel for
			//#endif
			for (long i = 0; i < sizeCells; i++) {
				CellHandle& cell = Tes.cellHandles[i];
				for (int j = 0; j < 4; j++) {
					//	VertexHandle vh = cell->vertex(j);
					if (cell->vertex(j)->info().isAlpha == true) {
						cell->info().p() = alphaBoundValue;
						cell->info().Pcondition = true;
						cell->info().isAlpha = true;
						alphaBoundingCells.push_back(cell);
					}
				}
			}
		}


		if (ppval && pxpos) applyUserDefinedPressure(Tri, *pxpos, *ppval);

		IPCells.clear();
		for (unsigned int n = 0; n < imposedP.size(); n++) {
			CellHandle cell = Tri.locate(CGT::Sphere(imposedP[n].first, 0));
			//check redundancy
			for (unsigned int kk = 0; kk < IPCells.size(); kk++) {
				if (cell == IPCells[kk]) cerr << "Two imposed pressures fall in the same cell." << endl;
				else if (cell->info().Pcondition)
					cerr << "Imposed pressure fall in a boundary condition." << endl;
			}
			IPCells.push_back(cell);
			cell->info().p() = imposedP[n].second;
			cell->info().Pcondition = true;
		}
		pressureChanged = false;

		IFCells.clear();
		for (unsigned int n = 0; n < imposedF.size(); n++) {
			CellHandle cell = Tri.locate(CGT::Sphere(imposedF[n].first, 0));
			//check redundancy
			for (unsigned int kk = 0; kk < IPCells.size(); kk++) {
				if (cell == IPCells[kk]) cerr << "Both flux and pressure are imposed in the same cell." << endl;
				else if (cell->info().Pcondition)
					cerr << "Imposed flux fall in a pressure boundary condition." << endl;
			}
			IFCells.push_back(cell);
			cell->info().Pcondition = false;
		}


		cavityCells.clear();
		for (unsigned int n = 0; n < imposedCavity.size(); n++) {
			CellHandle cell = Tri.locate(CGT::Sphere(imposedCavity[n], 0));
			cavityCells.push_back(cell);
			cell->info().isCavity = true;
			//cout << "cell set as cavity " << cell->info().id << endl;
		}
	}


	template <class Tesselation> void FlowBoundingSphere<Tesselation>::initializeTemperatures(Real tZero)
	{
		RTriangulation&     Tri = T[currentTes].Triangulation();
		FiniteCellsIterator cellEnd = Tri.finite_cells_end();

		for (FiniteCellsIterator cell = Tri.finite_cells_begin(); cell != cellEnd; cell++) {
			if (!cell->info().Tcondition && !cell->info().isGhost && !cell->info().blocked) cell->info().temp() = tZero;
		}
		for (int bound = 0; bound < 6; bound++) {
			int& id = *boundsIds[bound];
			thermalBoundingCells[bound].clear();
			if (id < 0) continue;
			ThermalBoundary& bi = thermalBoundary(id);
			if (!bi.fluxCondition) {
				VectorCell tmpCells;
				tmpCells.resize(10000);
				VCellIterator cells_it = tmpCells.begin();
				VCellIterator cells_end = Tri.incident_cells(T[currentTes].vertexHandles[id], cells_it);
				for (VCellIterator it = tmpCells.begin(); it != cells_end; it++) {
					(*it)->info().temp() = bi.value;
					(*it)->info().Tcondition = true;
					thermalBoundingCells[bound].push_back(*it);
				}
			}
		}
		//       if (ppval && pxpos) applyUserDefinedPressure(Tri,*pxpos,*ppval);

		//       ITCells.clear();
		//       for (unsigned int n=0; n<imposedT.size();n++) {
		//		CellHandle cell=Tri.locate(CGT::Sphere(imposedT[n].first,0));
		//check redundancy
		//		for (unsigned int kk=0;kk<IPCells.size();kk++){
		//			if (cell==IPCells[kk]) cerr<<"Two imposed pressures fall in the same cell."<<endl;
		//			else if  (cell->info().Pcondition) cerr<<"Imposed pressure fall in a boundary condition."<<endl;}
		//		ITCells.push_back(cell);
		//		cell->info().temp()=imposedT[n].second;
		//		cell->info().Tcondition=true;}
		//	pressureChanged=false;

		//	IFCells.clear();
		//	for (unsigned int n=0; n<imposedF.size();n++) {
		//		CellHandle cell=Tri.locate(CGT::Sphere(imposedF[n].first,0));
		//check redundancy
		//		for (unsigned int kk=0;kk<IPCells.size();kk++){
		//			if (cell==IPCells[kk]) cerr<<"Both flux and pressure are imposed in the same cell."<<endl;
		//			else if  (cell->info().Pcondition) cerr<<"Imposed flux fall in a pressure boundary condition."<<endl;}
		//		IFCells.push_back(cell);
		//		cell->info().Pcondition=false;}
	}


	template <class Tesselation> bool FlowBoundingSphere<Tesselation>::reApplyBoundaryConditions()
	{
		if (!pressureChanged) return false;
		for (int bound = 0; bound < 6; bound++) {
			int& id = *boundsIds[bound];
			if (id < 0) continue;
			Boundary& bi = boundary(id);
			if (!bi.flowCondition) {
				for (VCellIterator it = boundingCells[bound].begin(); it != boundingCells[bound].end(); it++) {
					(*it)->info().p() = bi.value;
					(*it)->info().Pcondition = true;
				}
			}
		}
		if (ppval && pxpos) applyUserDefinedPressure(T[currentTes].Triangulation(), *pxpos, *ppval);
		for (unsigned int n = 0; n < imposedP.size(); n++) {
			IPCells[n]->info().p() = imposedP[n].second;
			IPCells[n]->info().Pcondition = true;
		}
		pressureChanged = false;
		return true;
	}

	template <class Tesselation> void FlowBoundingSphere<Tesselation>::gaussSeidel(Real dt)
	{
		using math::max;
		using math::min;

		reApplyBoundaryConditions();
		RTriangulation& Tri = T[currentTes].Triangulation();
		int             j = 0;
		Real            m, n, dp_max, p_max, sum_p, p_moy, dp, sum_dp;
		Real            compFlowFactor = 0;
		vector<Real>    previousP;
		previousP.resize(Tri.number_of_finite_cells());
		const int num_threads = 1;
		bool      compressible = (fluidBulkModulus > 0);
#ifdef GS_OPEN_MP
		omp_set_num_threads(num_threads);
#endif

		if (debugOut) {
			cout << "tolerance = " << tolerance << endl;
			cout << "relax = " << relax << endl;
		}
		vector<Real> t_sum_p, t_dp_max, t_sum_dp, t_p_max;
		t_sum_dp.resize(num_threads);
		t_dp_max.resize(num_threads);
		t_p_max.resize(num_threads);
		t_sum_p.resize(num_threads);
		FiniteCellsIterator cellEnd = Tri.finite_cells_end();
#ifdef GS_OPEN_MP
		vector<FiniteCellsIterator> cells_its;
		for (FiniteCellsIterator cell = Tri.finite_cells_begin(); cell != cellEnd; cell++)
			if (!cell->info().Pcondition) cells_its.push_back(cell);
		int numCells = cells_its.size();
		cout << "cells_its.size() " << cells_its.size();
#endif
		// 	#pragma omp parallel shared(t_sum_dp, t_dp_max, sum_p, sum_dp,cells_its, j, Tri, relax)
		{
			do {
				int cell2 = 0;
				dp_max = 0;
				p_max = 0;
				p_moy = 0;
				sum_p = 0;
				sum_dp = 0;
#ifdef GS_OPEN_MP
				cell2 = numCells;
				for (int ii = 0; ii < num_threads; ii++)
					t_p_max[ii] = 0;
				for (int ii = 0; ii < num_threads; ii++)
					t_dp_max[ii] = 0;
				for (int ii = 0; ii < num_threads; ii++)
					t_sum_p[ii] = 0;
				for (int ii = 0; ii < num_threads; ii++)
					t_sum_dp[ii] = 0;
				int       kk = 0;
				const int numCells2 = numCells;
#pragma omp parallel for private(dp, m, n, kk) shared(tolerance, t_sum_dp, t_dp_max, sum_p, sum_dp, cells_its, j, Tri, relax) schedule(dynamic, 1000)
				for (kk = 0; kk < numCells2; kk++) {
					const FiniteCellsIterator& cell = cells_its[kk];
					{
#else
				int bb = -1;
				for (FiniteCellsIterator cell = Tri.finite_cells_begin(); cell != cellEnd; cell++) {
					bb++;
					if (!cell->info().Pcondition && !cell->info().blocked) {
						cell2++;
#endif
						if (compressible && j == 0) { previousP[bb] = cell->info().p(); }
						m = 0, n = 0;
						for (int j2 = 0; j2 < 4; j2++) {
							if (!Tri.is_infinite(cell->neighbor(j2))) {
								/// COMPRESSIBLE:
								if (compressible) {
									compFlowFactor = fluidBulkModulus * dt * cell->info().invVoidVolume();
									m += compFlowFactor * (cell->info().kNorm())[j2] * cell->neighbor(j2)->info().p();
									if (j == 0) n += compFlowFactor * (cell->info().kNorm())[j2];
								} else {
									/// INCOMPRESSIBLE
									m += (cell->info().kNorm())[j2] * cell->neighbor(j2)->info().p();
									if (math::isinf(m) && j < 10)
										cout << "(cell->info().kNorm())[j2] = " << (cell->info().kNorm())[j2]
										     << " cell->neighbor(j2)->info().p() = " << cell->neighbor(j2)->info().p()
										     << endl;
									if (j == 0) n += (cell->info().kNorm())[j2];
								}
							}
						}
						dp = cell->info().p();
						if (n != 0 || j != 0) {
							if (j == 0) {
								if (compressible) cell->info().invSumK = 1 / (1 + n);
								else
									cell->info().invSumK = 1 / n;
							}
							if (compressible) {
								/// COMPRESSIBLE cell->info().p() = ( (previousP - compFlowFactor*cell->info().dv()) + m ) / n ;
								cell->info().p()
								        = (((previousP[bb]
								             - ((fluidBulkModulus * dt * cell->info().invVoidVolume()) * (cell->info().dv())))
								            + m) * cell->info().invSumK
								           - cell->info().p())
								                * relax
								        + cell->info().p();
							} else {
								/// INCOMPRESSIBLE cell->info().p() =   - ( cell->info().dv() - m ) / ( n ) = ( -cell.info().dv() + m ) / n ;
								cell->info().p() = (-(cell->info().dv() - m) * cell->info().invSumK - cell->info().p()) * relax
								        + cell->info().p();
							}
#ifdef GS_OPEN_MP
#endif
						}
						dp -= cell->info().p();
#ifdef GS_OPEN_MP
						const int tn = omp_get_thread_num();
						t_sum_dp[tn] += math::abs(dp);
						t_dp_max[tn] = max(t_dp_max[tn], math::abs(dp));
						t_p_max[tn] = max(t_p_max[tn], math::abs(cell->info().p()));
						t_sum_p[tn] += math::abs(cell->info().p());
#else
						dp_max = max(dp_max, math::abs(dp));
						p_max = max(p_max, math::abs(cell->info().p()));
						sum_p += math::abs(cell->info().p());
						sum_dp += math::abs(dp);
#endif
					}
				}
#ifdef GS_OPEN_MP

				for (int ii = 0; ii < num_threads; ii++)
					p_max = max(p_max, t_p_max[ii]);
				for (int ii = 0; ii < num_threads; ii++)
					dp_max = max(dp_max, t_dp_max[ii]);
				for (int ii = 0; ii < num_threads; ii++)
					sum_p += t_sum_p[ii];
				for (int ii = 0; ii < num_threads; ii++)
					sum_dp += t_sum_dp[ii];
#endif
				p_moy = sum_p / cell2;

#ifdef GS_OPEN_MP
#pragma omp master
#endif
				j++;
#ifdef GS_OPEN_MP
			} while (j < 1500);
#else
			} while ((dp_max / p_max) > tolerance /*&& j<4000*/ /*&& ( dp_max > tolerance )*/ /* &&*/ /*( j<50 )*/);
#endif
		}
		if (debugOut) {
			cout << "pmax " << p_max << "; pmoy : " << p_moy << endl;
			cout << "iteration " << j << "; erreur : " << dp_max / p_max << endl;
		}
		computedOnce = true;
	}

	template <class Tesselation> Real FlowBoundingSphere<Tesselation>::boundaryFlux(unsigned int boundaryId)
	{
		if (noCache && T[!currentTes].Max_id() <= 0) return 0;
		bool            tes = noCache ? (!currentTes) : currentTes;
		RTriangulation& Tri = T[tes].Triangulation();
		Real            Q1 = 0;

		VectorCell tmpCells;
		tmpCells.resize(10000);
		VCellIterator cells_it = tmpCells.begin();

		VCellIterator cell_up_end = Tri.incident_cells(T[tes].vertexHandles[boundaryId], cells_it);
		for (VCellIterator it = tmpCells.begin(); it != cell_up_end; it++) {
			const CellHandle& cell = *it;
			if (cell->info().isGhost) continue;
			Q1 -= cell->info().dv();
			for (int j2 = 0; j2 < 4; j2++)
				Q1 += (cell->info().kNorm())[j2] * (cell->neighbor(j2)->info().shiftedP() - cell->info().shiftedP());
		}
		return Q1;
	}

	template <class Tesselation> Real FlowBoundingSphere<Tesselation>::boundaryArea(unsigned int boundaryId)
	{
		if (noCache && T[!currentTes].Max_id() <= 0) return 0;
		bool            tes = noCache ? (!currentTes) : currentTes;
		RTriangulation& Tri = T[tes].Triangulation();
		Real            A = 0;

		VectorCell tmpCells;
		tmpCells.resize(10000);
		VCellIterator cells_it = tmpCells.begin();

		VCellIterator cell_up_end = Tri.incident_cells(T[tes].vertexHandles[boundaryId], cells_it);
		for (VCellIterator it = tmpCells.begin(); it != cell_up_end; it++) {
			const CellHandle& cell = *it;
			if (cell->info().isGhost) continue;

			for (int j = 0; j < 4; j++) {
				if (cell->neighbor(j)->info().isFictious) continue;
				const CVector& Surfk = cell->info().facetSurfaces[j];
				Real           area = sqrt(Surfk.squared_length());
				A += cell->info().facetFluidSurfacesRatio[j] * area;
			}
		}
		return A;
	}

	template <class Tesselation> std::vector<std::vector<Real>> FlowBoundingSphere<Tesselation>::boundaryVel(unsigned int boundaryId)
	{
		//if (noCache && T[!currentTes].Max_id()<=0) return 0;
		bool            tes = noCache ? (!currentTes) : currentTes;
		RTriangulation& Tri = T[tes].Triangulation();
		//Real avgVel=0;
		//int numCells=0;
		std::vector<std::vector<Real>> velocities;
		VectorCell                     tmpCells;
		tmpCells.resize(10000);
		VCellIterator     cells_it = tmpCells.begin();
		std::vector<Real> vectorVel(4, 0);

		VCellIterator cell_up_end = Tri.incident_cells(T[tes].vertexHandles[boundaryId], cells_it);
		for (VCellIterator it = tmpCells.begin(); it != cell_up_end; it++) {
			//		numCells+=1;
			const CellHandle& cell = *it;
			if (cell->info().isGhost) continue;

			for (int j = 0; j < 4; j++) {
				if (cell->neighbor(j)->info().isFictious) continue;
				const CVector& Surfk = cell->info().facetSurfaces[j];
				Real           area = sqrt(Surfk.squared_length());
				vectorVel[0] = cell->info().averageVelocity()[0];
				vectorVel[1] = cell->info().averageVelocity()[1];
				vectorVel[2] = cell->info().averageVelocity()[2];
				vectorVel[3] = cell->info().facetFluidSurfacesRatio[j] * area;
				//vectorVel[4] = cell->info().kNorm()[j];
				velocities.push_back(vectorVel);
				//avgVel += sqrt(cell->info().averageVelocity().squared_length());
			}
		}
		//avgVel /= Real(numCells);
		return velocities;
	}


	template <class Tesselation> Real FlowBoundingSphere<Tesselation>::getCavityFlux()
	{
		Real cavityEdgeFlux = 0;
		//Real cavityVolume = 0;
		Tesselation& Tes = T[currentTes];
		const long   sizeCells = Tes.cellHandles.size();
#ifdef YADE_OPENMP
#pragma omp parallel for
#endif
		for (long i = 0; i < sizeCells; i++) {
			CellHandle& cell = Tes.cellHandles[i];
			if (!cell->info().isCavity || cell->info().isFictious || cell->info().blocked) continue; // only measuring cavity flux
			//cavityVolume += cell->info().volume();
			for (int j = 0; j < 4; j++) {
				if (cell->neighbor(j)->info().isCavity || cell->neighbor(j)->info().blocked) continue;
				cavityEdgeFlux
				        += -(cell->info().kNorm())[j] * (cell->info().p() - cell->neighbor(j)->info().p()); // influx to cavity will be negative
			}
		}

		return cavityEdgeFlux;
	}

	template <class Tesselation> void FlowBoundingSphere<Tesselation>::adjustCavityVolumeChange(Real dt, int stepsSinceLastMesh, Real /*pZero*/)
	{
		Real totalStep = dt * stepsSinceLastMesh;
		//Real Q1=0;
		netCavityFlux = 0;
		//Real cavityVolume = 0;
		Tesselation& Tes = T[currentTes];
		const long   sizeCells = Tes.cellHandles.size();
#ifdef YADE_OPENMP
#pragma omp parallel for
#endif
		for (long i = 0; i < sizeCells; i++) {
			CellHandle& cell = Tes.cellHandles[i];
			if (!cell->info().isCavity || cell->info().isFictious || cell->info().blocked) continue; // only measuring cavity flux
			//cavityVolume += cell->info().volume();
			for (int j = 0; j < 4; j++) {
				if (cell->neighbor(j)->info().isCavity || cell->neighbor(j)->info().blocked) continue;
				netCavityFlux
				        += -(cell->info().kNorm())[j] * (cell->info().p() - cell->neighbor(j)->info().p()); // influx to cavity will be negative
			}
		}

		// add imposed flux
		netCavityFlux += cavityFlux; // influx to cavity will be negative

		cavityDV += netCavityFlux / totalStep;
	}

	template <class Tesselation> void FlowBoundingSphere<Tesselation>::adjustCavityPressure(Real dt, int stepsSinceLastMesh, Real /*pZero*/)
	{
		Real totalStep = dt * stepsSinceLastMesh;
		//Real Q1=0;
		netCavityFlux = 0;
		Real         cavityVolume = 0;
		Tesselation& Tes = T[currentTes];
		const long   sizeCells = Tes.cellHandles.size();
#ifdef YADE_OPENMP
#pragma omp parallel for
#endif
		for (long i = 0; i < sizeCells; i++) {
			CellHandle& cell = Tes.cellHandles[i];
			if (!cell->info().isCavity || cell->info().isFictious || cell->info().blocked) continue; // only measuring cavity flux
			cavityVolume += cell->info().volume();
			for (int j = 0; j < 4; j++) {
				if (cell->neighbor(j)->info().isCavity || cell->neighbor(j)->info().blocked) continue;
				netCavityFlux += (cell->info().kNorm())[j]
				        * (cell->info().shiftedP() - cell->neighbor(j)->info().shiftedP()); // influx to cavity will be negative
			}
		}

		// add imposed flux
		netCavityFlux += cavityFlux; // influx to cavity will be negative

		// assign dynamic pcondition to cavity cells
		Real delP;
		if (cavityFluidDensity == 0) {
			delP = -netCavityFlux * totalStep / (equivalentCompressibility * cavityVolume);
		} // dv
		// or use density?
		else {
			Real cavityFluidDensityNew = ((cavityFluidDensity * cavityVolume) + (fluidRho * (-netCavityFlux * totalStep))) / (cavityVolume);
			delP = -(cavityFluidDensity / cavityFluidDensityNew - 1) * 1. / equivalentCompressibility;
			cavityFluidDensity = cavityFluidDensityNew;
		}
#ifdef YADE_OPENMP
#pragma omp parallel for
#endif
		for (long i = 0; i < sizeCells; i++) {
			CellHandle& cell = Tes.cellHandles[i];
			if (!cell->info().isCavity || cell->info().isFictious || cell->info().blocked) continue; // only measuring cavity flux
			cell->info().Pcondition = true;
			cell->info().p() += delP;
		}
		if (debugOut) cout << "flux added to cavity " << netCavityFlux << endl;
	}


	template <class Tesselation> void FlowBoundingSphere<Tesselation>::adjustCavityCompressibility(Real pZero)
	{
		Real sumP = 0;
		int  numCells = 0;
		netCavityFlux = 0;
		Tesselation& Tes = T[currentTes];
		const long   sizeCells = Tes.cellHandles.size();
#ifdef YADE_OPENMP
#pragma omp parallel for
#endif
		for (long i = 0; i < sizeCells; i++) {
			CellHandle& cell = Tes.cellHandles[i];
			if (!cell->info().isCavity || cell->info().isFictious || cell->info().blocked) continue;
			sumP += cell->info().p();
			numCells += 1;
			if (!controlCavityPressure) continue; // adjustCavPressure will track flux if active
			for (int j = 0; j < 4; j++) {
				if (cell->neighbor(j)->info().isCavity || cell->neighbor(j)->info().blocked) continue;
				netCavityFlux += (cell->info().kNorm())[j]
				        * (cell->info().shiftedP() - cell->neighbor(j)->info().shiftedP()); // reporting/debug purposes only, not used in model
			}
		}
		Real Pa = sumP / numCells;
		if (Pa == 0) {
			cerr << "0 pressure found while trying to account for air compressibility, invalid, setting to atmospheric" << endl;
			Pa = 1.0135e5;
		}
		Real Ca = 1. / Pa;
		Real phi = phiZero * (pZero / Pa);
		Real Cw = 1. / fluidBulkModulus;
		equivalentCompressibility = phi * Ca + (1. - phi) * Cw;
		if (debugOut) cout << "Equivalent compressibility " << equivalentCompressibility << endl;

		if (averageCavityPressure) {
#ifdef YADE_OPENMP
#pragma omp parallel for
#endif
			for (long i = 0; i < sizeCells; i++) {
				CellHandle& cell = Tes.cellHandles[i];
				if (!cell->info().isCavity || cell->info().isFictious || cell->info().blocked) continue;
				cell->info().p() = Pa; // set all cavity cell pressures the the average pressure...
			}
		}
	}


	template <class Tesselation> Real FlowBoundingSphere<Tesselation>::permeameter(Real PInf, Real PSup, Real Section, Real DeltaY, const char* file)
	{
		using math::max;
		using math::min;

		RTriangulation& Tri = T[currentTes].Triangulation();
		std::ofstream   kFile(file, std::ios::out);
		Real            Q2 = 0, Q1 = 0;
		int             cellQ1 = 0, cellQ2 = 0;
		Real            p_out_max = -10000000, p_out_min = 10000000, p_in_max = -100000000, p_in_min = 10000000, p_out_moy = 0, p_in_moy = 0;

		VectorCell tmpCells;
		tmpCells.resize(10000);
		VCellIterator cells_it = tmpCells.begin();

		VCellIterator cell_up_end = Tri.incident_cells(T[currentTes].vertexHandles[yMaxId], cells_it);
		for (VCellIterator it = tmpCells.begin(); it != cell_up_end; it++) {
			CellHandle& cell = *it;
			for (int j2 = 0; j2 < 4; j2++) {
				if (!cell->neighbor(j2)->info().Pcondition) {
					Q1 = Q1
					        + (cell->neighbor(j2)->info().kNorm())[Tri.mirror_index(cell, j2)]
					                * (cell->neighbor(j2)->info().p() - cell->info().p());
					cellQ1 += 1;
					p_out_max = max(cell->neighbor(j2)->info().p(), p_out_max);
					p_out_min = min(cell->neighbor(j2)->info().p(), p_out_min);
					p_out_moy += cell->neighbor(j2)->info().p();
				}
			}
		}

		VectorCell tmpCells2;
		tmpCells2.resize(10000);
		VCellIterator cells_it2 = tmpCells2.begin();

		VCellIterator cell_down_end = Tri.incident_cells(T[currentTes].vertexHandles[yMinId], cells_it2);
		for (VCellIterator it = tmpCells2.begin(); it != cell_down_end; it++) {
			CellHandle& cell = *it;
			for (int j2 = 0; j2 < 4; j2++) {
				if (!cell->neighbor(j2)->info().Pcondition) {
					Q2 = Q2
					        + (cell->neighbor(j2)->info().kNorm())[Tri.mirror_index(cell, j2)]
					                * (cell->info().p() - cell->neighbor(j2)->info().p());
					cellQ2 += 1;
					p_in_max = max(cell->neighbor(j2)->info().p(), p_in_max);
					p_in_min = min(cell->neighbor(j2)->info().p(), p_in_min);
					p_in_moy += cell->neighbor(j2)->info().p();
				}
			}
		}

		Real density = 1;
		//  declaration of ‘viscosity’ shadows a member of ‘yade::CGT::FlowBoundingSphere<_Tesselation>’ [-Werror=shadow]
		Real viscosity2 = viscosity;
		Real gravity = 1;
		Real Vdarcy = Q1 / Section;
		Real DeltaP = math::abs(PInf - PSup);
		Real DeltaH = DeltaP / (density * gravity);
		Real k = viscosity2 * Vdarcy * DeltaY / DeltaP; /**m²**/
		Real Ks = k * (density * gravity) / viscosity2; /**m/s**/

		if (debugOut) {
			cout << "the maximum superior pressure is = " << p_out_max << " the min is = " << p_out_min << endl;
			cout << "the maximum inferior pressure is = " << p_in_max << " the min is = " << p_in_min << endl;
			cout << "superior average pressure is " << p_out_moy / cellQ1 << endl;
			cout << "inferior average pressure is " << p_in_moy / cellQ2 << endl;
			cout << "celle comunicanti in basso = " << cellQ2 << endl;
			cout << "celle comunicanti in alto = " << cellQ1 << endl;
			cout << "The incoming average flow rate is = " << Q2 << " m^3/s " << endl;
			cout << "The outgoing average flow rate is = " << Q1 << " m^3/s " << endl;
			cout << "The gradient of charge is = " << DeltaH / DeltaY << " [-]" << endl;
			cout << "Darcy's velocity is = " << Vdarcy << " m/s" << endl;
			cout << "The permeability of the sample is = " << k << " m^2" << endl;
			cout << endl << "The hydraulic conductivity of the sample is = " << Ks << " m/s" << endl << endl;
		}
		kFile << "yMax id = " << yMaxId << "yMin id = " << yMinId << endl;
		kFile << "the maximum superior pressure is = " << p_out_max << " the min is = " << p_out_min << endl;
		kFile << "the maximum inferior pressure is = " << p_in_max << " the min is = " << p_in_min << endl;
		kFile << "superior average pressure is " << p_out_moy / cellQ2 << endl;
		kFile << "inferior average pressure is " << p_in_moy / cellQ1 << endl;
		kFile << "celle comunicanti in basso = " << cellQ2 << endl;
		kFile << "celle comunicanti in alto = " << cellQ1 << endl;
		kFile << "The incoming average flow rate is = " << Q2 << " m^3/s " << endl;
		kFile << "The outgoing average flow rate is = " << Q1 << " m^3/s " << endl;
		kFile << "The gradient of charge is = " << DeltaH / DeltaY << " [-]" << endl;
		kFile << "Darcy's velocity is = " << Vdarcy << " m/s" << endl;
		kFile << "The hydraulic conductivity of the sample is = " << Ks << " m/s" << endl;
		kFile << "The permeability of the sample is = " << k << " m^2" << endl;
		return k;
	}
	template <class Tesselation> void FlowBoundingSphere<Tesselation>::displayStatistics()
	{
		RTriangulation&     Tri = T[currentTes].Triangulation();
		int                 Zero = 0, Inside = 0, Fictious = 0;
		FiniteCellsIterator cellEnd = Tri.finite_cells_end();
		for (FiniteCellsIterator cell = Tri.finite_cells_begin(); cell != cellEnd; cell++) {
			int zeros = 0;
			for (int j = 0; j != 4; j++) {
				if ((cell->info().kNorm())[j] == 0) { zeros += 1; }
			}
			if (zeros == 4) { Zero += 1; }
			if (!cell->info().fictious()) {
				Inside += 1;
			} else {
				Fictious += 1;
			}
		}
		int fict = 0, real = 0;
		for (FiniteVerticesIterator v = Tri.finite_vertices_begin(); v != Tri.finite_vertices_end(); ++v) {
			if (v->info().isFictious) fict += 1;
			else
				real += 1;
		}
		long Vertices = Tri.number_of_vertices();
		long Cells = Tri.number_of_finite_cells();
		long Facets = Tri.number_of_finite_facets();
		if (debugOut) {
			cout << "zeros = " << Zero << endl;
			cout << "There are " << Vertices << " vertices, dont " << fict << " fictious et " << real << " reeeeeel" << std::endl;
			cout << "There are " << Cells << " cells " << std::endl;
			cout << "There are " << Facets << " facets " << std::endl;
			cout << "There are " << Inside << " cells INSIDE." << endl;
			cout << "There are " << Fictious << " cells FICTIOUS." << endl;
		}

		num_particles = real;
	}


	template <class Tesselation> void FlowBoundingSphere<Tesselation>::saveVtk(const char* folder, bool withBoundaries)
	{
		vector<int>
		            allIds; //an ordered list of cell ids (from begin() to end(), for vtk table lookup), some ids will appear multiple times since boundary cells are splitted into multiple tetrahedra
		vector<int> fictiousN;


		static unsigned int number = 0;
		char                filename[250];
		mkdir(folder, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		sprintf(filename, "%s/out_%d.vtk", folder, number++);

		basicVTKwritter vtkWrite(0, 0);
		saveMesh(vtkWrite, withBoundaries, allIds, fictiousN, filename);

		// Right after remeshing "currentTes" points to the new mesh which doesn't contain valid pressure yet, this is detected by noCache=true, in such case get data from previous mesh
		RTriangulation& Tri = T[noCache ? (!currentTes) : currentTes].Triangulation();
		VectorCell&     cellHandles = T[noCache ? (!currentTes) : currentTes].cellHandles;

		if (permeabilityMap) {
			vtkWrite.begin_data("Permeability", CELL_DATA, SCALARS, FLOAT);
			for (FiniteCellsIterator cell = Tri.finite_cells_begin(); cell != Tri.finite_cells_end(); ++cell) {
				bool isDrawable = cell->info().isReal() && cell->vertex(0)->info().isReal() && cell->vertex(1)->info().isReal()
				        && cell->vertex(2)->info().isReal() && cell->vertex(3)->info().isReal();
				if (isDrawable) {
					vtkWrite.write_data(
					        (cell->info().kNorm()[0] + cell->info().kNorm()[1] + cell->info().kNorm()[2] + cell->info().kNorm()[3]) / 4.);
				}
			}
			vtkWrite.end_data();
		} else { //normal case
			vtkWrite.begin_data("Pressure", CELL_DATA, SCALARS, FLOAT);
			for (unsigned kk = 0; kk < allIds.size(); kk++)
				vtkWrite.write_data(cellHandles[allIds[kk]]->info().p());
			vtkWrite.end_data();

			if (thermalEngine) {
				vtkWrite.begin_data("Temperature", CELL_DATA, SCALARS, FLOAT);
				for (unsigned kk = 0; kk < allIds.size(); kk++) {
					bool isDrawable = cellHandles[allIds[kk]]->info().isReal() && cellHandles[allIds[kk]]->vertex(0)->info().isReal()
					        && cellHandles[allIds[kk]]->vertex(1)->info().isReal() && cellHandles[allIds[kk]]->vertex(2)->info().isReal()
					        && cellHandles[allIds[kk]]->vertex(3)->info().isReal();
					if (isDrawable) vtkWrite.write_data(cellHandles[allIds[kk]]->info().temp());
				}
				vtkWrite.end_data();

				vtkWrite.begin_data("Reynolds", CELL_DATA, SCALARS, FLOAT);
				for (unsigned kk = 0; kk < allIds.size(); kk++) {
					bool isDrawable = cellHandles[allIds[kk]]->info().isReal() && cellHandles[allIds[kk]]->vertex(0)->info().isReal()
					        && cellHandles[allIds[kk]]->vertex(1)->info().isReal() && cellHandles[allIds[kk]]->vertex(2)->info().isReal()
					        && cellHandles[allIds[kk]]->vertex(3)->info().isReal();
					if (isDrawable) vtkWrite.write_data(cellHandles[allIds[kk]]->info().Reynolds);
				}
				vtkWrite.end_data();

				vtkWrite.begin_data("Tcondition", CELL_DATA, SCALARS, FLOAT);
				for (FiniteCellsIterator cell = Tri.finite_cells_begin(); cell != Tri.finite_cells_end(); ++cell) {
					bool isDrawable = cell->info().isReal() && cell->vertex(0)->info().isReal() && cell->vertex(1)->info().isReal()
					        && cell->vertex(2)->info().isReal() && cell->vertex(3)->info().isReal();
					if (isDrawable) { vtkWrite.write_data(cell->info().Tcondition); }
				}
				vtkWrite.end_data();
			}

			vtkWrite.begin_data("cavity", CELL_DATA, SCALARS, FLOAT);
			for (FiniteCellsIterator cell = Tri.finite_cells_begin(); cell != Tri.finite_cells_end(); ++cell) {
				bool isDrawable = cell->info().isReal() && cell->vertex(0)->info().isReal() && cell->vertex(1)->info().isReal()
				        && cell->vertex(2)->info().isReal() && cell->vertex(3)->info().isReal();
				if (isDrawable) { vtkWrite.write_data(cell->info().isCavity); }
			}
			vtkWrite.end_data();

			vtkWrite.begin_data("alpha", CELL_DATA, SCALARS, FLOAT);
			for (FiniteCellsIterator cell = Tri.finite_cells_begin(); cell != Tri.finite_cells_end(); ++cell) {
				bool isDrawable = cell->info().isReal() && cell->vertex(0)->info().isReal() && cell->vertex(1)->info().isReal()
				        && cell->vertex(2)->info().isReal() && cell->vertex(3)->info().isReal();
				if (isDrawable) { vtkWrite.write_data(cell->info().isAlpha); }
			}
			vtkWrite.end_data();

			vtkWrite.begin_data("Pcondition", CELL_DATA, SCALARS, FLOAT);
			for (FiniteCellsIterator cell = Tri.finite_cells_begin(); cell != Tri.finite_cells_end(); ++cell) {
				bool isDrawable = cell->info().isReal() && cell->vertex(0)->info().isReal() && cell->vertex(1)->info().isReal()
				        && cell->vertex(2)->info().isReal() && cell->vertex(3)->info().isReal();
				if (isDrawable) { vtkWrite.write_data(cell->info().Pcondition); }
			}
			vtkWrite.end_data();

			vtkWrite.begin_data("fictious", CELL_DATA, SCALARS, INT);
			for (unsigned kk = 0; kk < allIds.size(); kk++)
				vtkWrite.write_data(fictiousN[kk]);
			vtkWrite.end_data();

			vtkWrite.begin_data("id", CELL_DATA, SCALARS, INT);
			for (unsigned kk = 0; kk < allIds.size(); kk++)
				vtkWrite.write_data(allIds[kk]);
			vtkWrite.end_data();

			averageRelativeCellVelocity();
			vtkWrite.begin_data("Velocity", CELL_DATA, VECTORS, FLOAT);
			for (unsigned kk = 0; kk < allIds.size(); kk++)
				vtkWrite.write_data(
				        cellHandles[allIds[kk]]->info().averageVelocity()[0],
				        cellHandles[allIds[kk]]->info().averageVelocity()[1],
				        cellHandles[allIds[kk]]->info().averageVelocity()[2]);
			vtkWrite.end_data();
		}
		vtkWrite.close();
	}

	template <class Tesselation>
	void FlowBoundingSphere<Tesselation>::saveMesh(
	        basicVTKwritter& vtkfile, bool withBoundaries, vector<int>& allIds, vector<int>& fictiousN, const char* filename)
	{
		using math::max;
		using math::min;

		/*
	Most of this code is actually to handle the special cases along the bounding walls/spheres, where irregular polyhedra needs to be decomposed in many tetrahedra for vtk display.
	- fictious=1: take a triangular (pseudo-)prism formed by 3 spheres and 3 projections and split it in 3 tetrahedra
	- fictious=2: take a rectangular (pseudo-)prism formed by 2 spheres and 6 projections and split it in 6 tetrahedra
	- fictious=3: take a rectangular (pseudo-)prism formed by 1 sphere and 7 projections and split it in 6 tetrahedra
	"extraVertices" and "extraPoly" contains the additional points and polyhedra.

	*/

		if (noCache && T[!currentTes].Max_id() <= 0) {
			cout << "Triangulation does not exist. Sorry." << endl;
			return;
		}
		RTriangulation& Tri = T[noCache ? (!currentTes) : currentTes].Triangulation();
		int             firstReal = -1;

		vector<pair<vector<int>, unsigned>> extraPoly;     //for fictious cells on the boundaries
		vector<CGT::CVector>                extraVertices; //for projections of real vertices on the boundaries

		vector<int> fictiousNExtra;
		//count fictious vertices and cells
		vtkInfiniteVertices = vtkInfiniteCells = 0;
		unsigned extraV = 0, extraCells = 0;

		FiniteCellsIterator cellEnd = Tri.finite_cells_end();
		for (FiniteCellsIterator cell = Tri.finite_cells_begin(); cell != cellEnd; cell++) {
			bool fictV[4];
			for (unsigned k = 0; k < 4; k++)
				fictV[k] = !cell->vertex(k)->info().isReal();
			unsigned fict = unsigned(fictV[0]) + unsigned(fictV[1]) + unsigned(fictV[2]) + unsigned(fictV[3]);
			if (fict == 0) {
				allIds.push_back(cell->info().id);
				fictiousN.push_back(0);
			} else {
				vtkInfiniteCells += 1;
				if (withBoundaries) {
					int newId = -(extraVertices.size() + 1); //negative ids starting from -1 for fictious vertices
					switch (fict) {
						case 1: {
							int k = 0;
							for (; k < 4; k++)
								if (fictV[k]) break; //find the fictious one
							vector<int> poly;
							for (int j = 0; j < 4; j++) {
								if (!fictV[j]) {
									Real dist = (cell->vertex(k)->point().point() - cell->vertex(j)->point().point())
									                    .squared_length();
									dist = 1. - sqrt(cell->vertex(k)->point().weight() / dist);
									CVector C1C2
									        = dist * (cell->vertex(k)->point().point() - cell->vertex(j)->point().point());
									extraVertices.push_back((cell->vertex(j)->point().point() + C1C2) - CGAL::ORIGIN);
									poly.push_back(cell->vertex(j)->info().id());
								}
							}
							// then the real ones
							poly.push_back(newId--);
							poly.push_back(newId--);
							poly.push_back(newId--);
							// store the extra poly for later use and update counters
							extraPoly.push_back(pair<vector<int>, int>(poly, cell->info().id)); //
							for (int kk = 0; kk < 3; kk++)
								fictiousNExtra.push_back(1);
							extraCells += 3;
							extraV += 3;
						} break;
						case 2: {
							unsigned k = 0;
							unsigned k1 = 0, k2 = 0, s1 = 0, s2 = 0;
							for (; k < 4; k++)
								if (fictV[k]) {
									k1 = k;
									k++;
									break;
								} //find the first fictious
							for (; k < 4; k++)
								if (fictV[k]) {
									k2 = k;
									k++;
									break;
								} //find the second fictious
							for (k = 0; k < 4; k++)
								if (!fictV[k]) {
									s1 = k;
									k++;
									break;
								} //find the first real
							for (; k < 4; k++)
								if (!fictV[k]) {
									s2 = k;
									k++;
									break;
								} //find the second real

							CVector c1 = cell->vertex(s1)->point().point() - CGAL::ORIGIN;
							CVector c2 = cell->vertex(s2)->point().point() - CGAL::ORIGIN;
							CVector v11 = cell->vertex(k1)->point().point() - cell->vertex(s1)->point().point();
							CVector v21 = cell->vertex(k1)->point().point() - cell->vertex(s2)->point().point();
							CVector v12 = cell->vertex(k2)->point().point() - cell->vertex(s1)->point().point();
							CVector v22 = cell->vertex(k2)->point().point() - cell->vertex(s2)->point().point();
							CVector v13;
							CVector v23;

							Real dist11 = v11.squared_length();
							dist11 = 1. - sqrt(cell->vertex(k1)->point().weight() / dist11);
							v11 = v11 * dist11;
							Real dist21 = v21.squared_length();
							dist21 = 1. - sqrt(cell->vertex(k1)->point().weight() / dist21);
							v21 = v21 * dist21;
							Real dist12 = v12.squared_length();
							dist12 = 1. - sqrt(cell->vertex(k2)->point().weight() / dist12);
							v12 = v12 * dist12;
							Real dist22 = v22.squared_length();
							dist22 = 1. - sqrt(cell->vertex(k2)->point().weight() / dist22);
							v22 = v22 * dist22;
							v13 = v12 + v11;
							v23 = v22 + v21;
							vector<int> poly1, poly2; //split the cube in two prisms
							poly1.push_back(cell->vertex(s1)->info().id());
							poly2.push_back(cell->vertex(s1)->info().id());
							extraVertices.push_back(c1 + v13);
							poly1.push_back(newId);
							poly2.push_back(newId--);
							extraVertices.push_back(c1 + v11);
							poly1.push_back(newId--);
							extraVertices.push_back(c1 + v12);
							poly2.push_back(newId--);
							poly1.push_back(cell->vertex(s2)->info().id());
							poly2.push_back(cell->vertex(s2)->info().id());
							extraVertices.push_back(c2 + v23);
							poly1.push_back(newId);
							poly2.push_back(newId--);
							extraVertices.push_back(c2 + v21);
							poly1.push_back(newId--);
							extraVertices.push_back(c2 + v22);
							poly2.push_back(newId--);
							extraPoly.push_back(pair<vector<int>, int>(poly1, cell->info().id));
							extraPoly.push_back(pair<vector<int>, int>(poly2, cell->info().id));
							for (int kk = 0; kk < 6; kk++)
								fictiousNExtra.push_back(2);
							extraCells += 6;
							extraV += 6;
						} break;

						case 3: {
							unsigned s1 = 0, k1 = 0, k2 = 0, k3 = 0;
							for (unsigned k = 0; k < 4; k++)
								if (!fictV[k]) {
									s1 = k;
									k++;
									break;
								}                  //find the first real
							k1 = facetVertices[s1][0]; //handy...
							k3 = facetVertices[s1][1];
							k2 = facetVertices[s1][2];

							CVector c = cell->vertex(s1)->point().point() - CGAL::ORIGIN;
							CVector v11 = cell->vertex(k1)->point().point() - cell->vertex(s1)->point().point();
							CVector v12 = cell->vertex(k2)->point().point() - cell->vertex(s1)->point().point();
							CVector v13 = cell->vertex(k3)->point().point() - cell->vertex(s1)->point().point();

							Real dist11 = v11.squared_length();
							dist11 = 1. - sqrt(cell->vertex(k1)->point().weight() / dist11);
							v11 = v11 * dist11;
							Real dist12 = v12.squared_length();
							dist12 = 1. - sqrt(cell->vertex(k2)->point().weight() / dist12);
							v12 = v12 * dist12;
							Real dist13 = v13.squared_length();
							dist13 = 1. - sqrt(cell->vertex(k3)->point().weight() / dist13);
							v13 = v13 * dist13;
							CVector v1v2 = v11 + v12;

							vector<int> poly1, poly2; //split the cube in two prisms
							poly1.push_back(cell->vertex(s1)->info().id());
							poly2.push_back(cell->vertex(s1)->info().id());
							extraVertices.push_back(c + v1v2);
							poly1.push_back(newId);
							poly2.push_back(newId--);
							extraVertices.push_back(c + v11);
							poly1.push_back(newId--);
							extraVertices.push_back(c + v12);
							poly2.push_back(newId--);
							extraVertices.push_back(c + v13);
							poly1.push_back(newId);
							poly2.push_back(newId--);
							extraVertices.push_back(c + v1v2 + v13);
							poly1.push_back(newId);
							poly2.push_back(newId--);
							extraVertices.push_back(c + v11 + v13);
							poly1.push_back(newId--);
							extraVertices.push_back(c + v12 + v13);
							poly2.push_back(newId--);

							extraPoly.push_back(pair<vector<int>, int>(poly1, cell->info().id));
							extraPoly.push_back(pair<vector<int>, int>(poly2, cell->info().id));
							for (int kk = 0; kk < 6; kk++)
								fictiousNExtra.push_back(3);
							extraCells += 6;
							extraV += 7;
						} break;
						default: throw std::runtime_error(__FILE__ " : switch default case error.");
					}
				}
			}
		}
		fictiousN.insert(fictiousN.end(), fictiousNExtra.begin(), fictiousNExtra.end());

		for (FiniteVerticesIterator v = Tri.finite_vertices_begin(); v != Tri.finite_vertices_end(); ++v) {
			if (!v->info().isReal()) vtkInfiniteVertices += 1;
			else if (firstReal == -1)
				firstReal = vtkInfiniteVertices;
		}

		vtkfile.setNums(
		        (unsigned int)Tri.number_of_vertices() - vtkInfiniteVertices + extraV,
		        (unsigned int)Tri.number_of_finite_cells() - vtkInfiniteCells + extraCells);
		vtkfile.open(filename, "Generated by Yade-dem.org");

		//!TEMPORARY FIX:
		//paraview needs zero-based vertex ids (from 0 ... numRealVertices)
		//in presence of clumps vertex ids are not zero-based anymore
		//to fix the vkt output vertex ids will be replaced by zero-based ones (CAUTION: output vertex ids != Yade vertex ids!)

		std::map<int, int> vertexIdMap;
		int                numVertices = 0;
		unsigned int       minId = 1;

		vtkfile.begin_vertices();
		Real x, y, z;
		for (FiniteVerticesIterator v = Tri.finite_vertices_begin(); v != Tri.finite_vertices_end(); ++v) {
			if (v->info().isReal()) {
				x = (Real)(v->point().point()[0]);
				y = (Real)(v->point().point()[1]);
				z = (Real)(v->point().point()[2]);
				vtkfile.write_point(x, y, z);
				vertexIdMap[v->info().id() - firstReal] = numVertices;
				minId = min(minId, v->info().id() - firstReal);
				numVertices += 1;
			}
		}
		// now the extra vertices
		int extra = -1;
		int numAllVertices = numVertices;
		for (auto pt = extraVertices.begin(); pt != extraVertices.end(); pt++) {
			vtkfile.write_point((Real)(*pt)[0], (Real)(*pt)[1], (Real)(*pt)[2]);
			vertexIdMap[extra--] = numAllVertices++;
		}
		vtkfile.end_vertices();

		vtkfile.begin_cells();
		int id0, id1, id2, id3;
		for (FiniteCellsIterator cell = Tri.finite_cells_begin(); cell != Tri.finite_cells_end(); ++cell) {
			bool isDrawable = cell->info().isReal() && cell->vertex(0)->info().isReal() && cell->vertex(1)->info().isReal()
			        && cell->vertex(2)->info().isReal() && cell->vertex(3)->info().isReal();
			if (isDrawable) {
				id0 = cell->vertex(0)->info().id() - firstReal;
				id1 = cell->vertex(1)->info().id() - firstReal;
				id2 = cell->vertex(2)->info().id() - firstReal;
				id3 = cell->vertex(3)->info().id() - firstReal;
				if (minId != 0) vtkfile.write_cell(vertexIdMap[id0], vertexIdMap[id1], vertexIdMap[id2], vertexIdMap[id3]);
				else
					vtkfile.write_cell(id0, id1, id2, id3);
			}
		}
		// now the extra tetrahedra
		int id4, id5;
		for (auto p = extraPoly.begin(); p != extraPoly.end(); p++) {
			switch (p->first.size()) {
				case 6:
					allIds.push_back(p->second);
					allIds.push_back(p->second);
					allIds.push_back(p->second);
					if (minId != 0) {
						id0 = p->first[0];
						id1 = p->first[1];
						id2 = p->first[2];
						id3 = p->first[3];
						id4 = p->first[4];
						id5 = p->first[5]; // extra vertices, negative ids
						vtkfile.write_cell(vertexIdMap[id0], vertexIdMap[id3], vertexIdMap[id4], vertexIdMap[id5]);
						vtkfile.write_cell(vertexIdMap[id0], vertexIdMap[id5], vertexIdMap[id1], vertexIdMap[id2]);
						vtkfile.write_cell(vertexIdMap[id0], vertexIdMap[id1], vertexIdMap[id4], vertexIdMap[id5]);
					} else {
						//real spheres ave positive ids, fictious ones ahave negative, that's how we guess the rank in the liste
						id0 = p->first[0] > 0 ? p->first[0] - firstReal : numVertices - p->first[0] - 1;
						id1 = p->first[1] > 0 ? p->first[1] - firstReal : numVertices - p->first[1] - 1;
						id2 = p->first[2] > 0 ? p->first[2] - firstReal : numVertices - p->first[2] - 1;
						id3 = p->first[3] > 0 ? p->first[3] - firstReal : numVertices - p->first[3] - 1;
						id4 = p->first[4] > 0 ? p->first[4] - firstReal : numVertices - p->first[4] - 1;
						id5 = p->first[5] > 0 ? p->first[5] - firstReal : numVertices - p->first[5] - 1;
						vtkfile.write_cell(id0, id3, id4, id5);
						vtkfile.write_cell(id0, id5, id1, id2);
						vtkfile.write_cell(id0, id1, id4, id5);
					}
					break;
				default: throw std::runtime_error(__FILE__ " : switch default case error.");
			}
		}
		vtkfile.end_cells();
	}

#ifdef XVIEW
	template <class Tesselation> void FlowBoundingSphere<Tesselation>::dessineTriangulation(Vue3D& Vue, RTriangulation& T)
	{
		Real* Segments = NULL;
		long  N_seg = newListeEdges(T, &Segments);
		Vue.Dessine_Segment(Segments, N_seg);
		deleteListeEdges(&Segments, N_seg);
	}
	template <class Tesselation> void FlowBoundingSphere<Tesselation>::dessineShortTesselation(Vue3D& Vue, Tesselation& Tes)
	{
		if (!Tes.computed()) Tes.compute();
		Real* Segments = NULL;
		long  N_seg = Tes.newListeShortEdges(&Segments);
		Vue.Dessine_Segment(Segments, N_seg);
		deleteListeEdges(&Segments, N_seg);
	}
#endif
	template <class Tesselation> void FlowBoundingSphere<Tesselation>::generateVoxelFile()
	{
		RTriangulation& Tri = T[currentTes].Triangulation();
		Real            l = 1;
		int             dx = 200;
		Real            eps = l / dx;

		std::ofstream voxelfile("MATRIX", std::ios::out);
		bool          solid = false;

		for (Real y = 0; y <= l; y += eps) {
			for (Real z = 0; z <= l; z += eps) {
				for (Real x = 0; x <= l; x += eps) {
					solid = false;

					for (FiniteVerticesIterator vIt = Tri.finite_vertices_begin(); vIt != Tri.finite_vertices_end(); vIt++) {
						Real radius = sqrt(vIt->point().weight());
						if ((sqrt(pow((x - (vIt->point()[0])), 2) + pow((y - (vIt->point()[1])), 2) + pow((z - (vIt->point()[2])), 2)))
						    <= radius)
							solid = true;
					}
					if (solid) voxelfile << 1;
					else
						voxelfile << 0;
				}
				voxelfile << endl;
			}
		}
	}

	template <class Tesselation> void FlowBoundingSphere<Tesselation>::printVertices()
	{
		RTriangulation& Tri = T[currentTes].Triangulation();
		std::ofstream   file("vertices.txt", std::ios::out);
		file << "id x y z r alpha fictious" << endl;
		for (FiniteVerticesIterator vIt = Tri.finite_vertices_begin(); vIt != Tri.finite_vertices_end(); vIt++) {
			file << vIt->info().id() << " " << vIt->point()[0] << " " << vIt->point()[1] << " " << vIt->point()[2] << " "
			     << " " << sqrt(vIt->point().weight()) << " " << vIt->info().isAlpha << " " << vIt->info().isFictious << endl;
		}
		file.close();
	}

	template <class Tesselation>
	Real FlowBoundingSphere<Tesselation>::samplePermeability(Real& xMin2, Real& xMax2, Real& yMin2, Real& yMax2, Real& zMin2, Real& zMax2 /*, string key*/)
	{
		Real Section = (xMax2 - xMin2) * (zMax2 - zMin2);
		Real DeltaY = yMax2 - yMin2;
		boundary(yMinId).flowCondition = 0;
		boundary(yMaxId).flowCondition = 0;
		boundary(yMinId).value = 0;
		boundary(yMaxId).value = 1;
		Real pZero = math::abs((boundary(yMinId).value - boundary(yMaxId).value) / 2);
		initializePressure(pZero);
		gaussSeidel();
		const char* kk = "Permeability";
		return permeameter(boundary(yMinId).value, boundary(yMaxId).value, Section, DeltaY, kk);
	}
	template <class Tesselation> bool FlowBoundingSphere<Tesselation>::isInsideSphere(Real& x, Real& y, Real& z)
	{
		RTriangulation& Tri = T[currentTes].Triangulation();
		for (FiniteVerticesIterator vIt = Tri.finite_vertices_begin(); vIt != Tri.finite_vertices_end(); vIt++) {
			Real radius = vIt->point().weight();
			if (pow((x - (vIt->point()[0])), 2) + pow((y - (vIt->point()[1])), 2) + pow((z - (vIt->point()[2])), 2) <= radius) return true;
		}
		return false;
	}
	template <class Tesselation> void FlowBoundingSphere<Tesselation>::sliceField(const char* filename)
	{
		/** Pressure field along one cutting plane **/
		RTriangulation& Tri = T[noCache ? (!currentTes) : currentTes].Triangulation();
		CellHandle      permeameter;

		std::ofstream consFile(filename, std::ios::out);

		int  intervals = 400;
		Real Ry = (yMax - yMin) / intervals;
		Real Rz = (zMax - zMin) / intervals;
		Real X = 0.5;
		for (Real Y = min(yMax, yMin); Y <= max(yMax, yMin); Y = Y + math::abs(Ry)) {
			for (Real Z = min(zMin, zMax); Z <= max(zMin, zMax); Z = Z + math::abs(Rz)) {
				permeameter = Tri.locate(Point(X, Y, Z));
				consFile << permeameter->info().p() << " ";
			}
			consFile << endl;
		}
		consFile << endl;
		consFile.close();
	}

	template <class Tesselation> void FlowBoundingSphere<Tesselation>::computeEdgesSurfaces()
	{
		using math::max;
		using math::min;

		RTriangulation& Tri = T[currentTes].Triangulation();
		//first, copy interacting pairs and normal lub forces form prev. triangulation in a sorted structure for initializing the new lub. Forces
		vector<vector<pair<unsigned int, Real>>> lubPairs;
		lubPairs.resize(Tri.number_of_vertices() + 1);
		for (unsigned int k = 0; k < edgeNormalLubF.size(); k++)
			lubPairs[min(edgeIds[k].first->id(), edgeIds[k].second->id())].push_back(
			        pair<int, Real>(max(edgeIds[k].first->id(), edgeIds[k].second->id()), edgeNormalLubF[k]));

		//Now we reset the containers and initialize them
		edgeSurfaces.clear();
		edgeIds.clear();
		edgeNormalLubF.clear();
		//	FiniteEdgesIterator ed_it;
		for (FiniteEdgesIterator ed_it = Tri.finite_edges_begin(); ed_it != Tri.finite_edges_end(); ed_it++) {
			const VertexInfo& vi1 = (ed_it->first)->vertex(ed_it->second)->info();
			const VertexInfo& vi2 = (ed_it->first)->vertex(ed_it->third)->info();

			//We eliminate edges that would be periodic replications or involving two bounding objects, i.e. the min id must be non-ghost and non-fictious
			if (vi1.id() < vi2.id()) {
				if (vi1.isFictious || vi2.isGhost) continue;
			} else if (vi2.isFictious || vi2.isGhost)
				continue;
			Real area = T[currentTes].computeVFacetArea(ed_it);
			edgeSurfaces.push_back(area);
			unsigned int id1 = vi1.id();
			unsigned int id2 = vi2.id();
			edgeIds.push_back(pair<const VertexInfo*, const VertexInfo*>(&vi1, &vi2));
			//For persistant edges, we must transfer the lub. force value from the older triangulation structure
			if (id1 > id2) swap(id1, id2);
			unsigned int i = 0;
			//Look for the pair (id1,id2) in lubPairs
			while (i < lubPairs[id1].size()) {
				if (lubPairs[id1][i].first == id2) {
					//it's found, we copy the lub force
					edgeNormalLubF.push_back(lubPairs[id1][i].second);
					break;
				}
				++i;
			}
			// not found, we initialize with zero lub force
			if (i == lubPairs[id1].size()) edgeNormalLubF.push_back(0);
		}
	}

	template <class Tesselation>
	Vector3r FlowBoundingSphere<Tesselation>::computeViscousShearForce(const Vector3r& deltaV, const int& edge_id, const Real& Rh)
	{
		Vector3r tau = deltaV * viscosity / Rh;
		return tau * edgeSurfaces[edge_id];
	}

	template <class Tesselation>
	Vector3r FlowBoundingSphere<Tesselation>::computeShearLubricationForce(
	        const Vector3r& deltaShearV, const Real& dist, const int& /*edge_id*/, const Real& eps, const Real& centerDist, const Real& meanRad)
	{
		Real     d = math::max(dist, 0.) + 2. * eps * meanRad;
		Vector3r viscLubF = 0.5 * Mathr::PI * viscosity * (-2 * meanRad + centerDist * log(centerDist / d)) * deltaShearV;
		return viscLubF;
	}

	template <class Tesselation>
	Vector3r FlowBoundingSphere<Tesselation>::computePumpTorque(
	        const Vector3r& deltaShearAngV, const Real& dist, const int& /*edge_id*/, const Real& eps, const Real& meanRad)
	{
		Real     d = math::max(dist, 0.) + 2. * eps * meanRad;
		Vector3r viscPumpC = Mathr::PI * viscosity * pow(meanRad, 3) * (3. / 20. * log(meanRad / d) + 63. / 500. * (d / meanRad) * log(meanRad / d))
		        * deltaShearAngV;
		return viscPumpC;
	}

	template <class Tesselation>
	Vector3r FlowBoundingSphere<Tesselation>::computeTwistTorque(
	        const Vector3r& deltaNormAngV, const Real& dist, const int& /*edge_id*/, const Real& eps, const Real& meanRad)
	{
		Real     d = math::max(dist, 0.) + 2. * eps * meanRad;
		Vector3r twistC = Mathr::PI * viscosity * pow(meanRad, 2) * d * log(meanRad / d) * deltaNormAngV;
		return twistC;
	}


	template <class Tesselation>
	Real FlowBoundingSphere<Tesselation>::computeNormalLubricationForce(
	        const Real& deltaNormV, const Real& dist, const int& edge_id, const Real& eps, const Real& stiffness, const Real& dt, const Real& meanRad)
	{
		//FIXME: here introduce elasticity
		Real d = math::max(dist, 0.) + 2. * eps * meanRad; //account for grains roughness
		if (stiffness > 0) {
			const Real k = stiffness * meanRad;
			Real       prevForce = edgeNormalLubF[edge_id];
			Real       instantVisc = 1.5 * Mathr::PI * pow(meanRad, 2) * viscosity / (d - prevForce / k);
			Real       normLubF = instantVisc * (deltaNormV + prevForce / (k * dt)) / (1 + instantVisc / (k * dt));
			edgeNormalLubF[edge_id] = normLubF;
			return normLubF;
		} else {
			Real normLubF = (1.5 * Mathr::PI * pow(meanRad, 2) * viscosity * deltaNormV) / d;
			return normLubF;
		}
	}

	template <class Tesselation> Real FlowBoundingSphere<Tesselation>::fractionalSolidArea(CellHandle cell, int j)
	{
		Real area;
		int  k = 0, l = 0, m = 0;
		if (j == 0) {
			k = 1;
			l = 2;
			m = 3;
		}
		if (j == 1) {
			k = 0;
			l = 2;
			m = 3;
		}
		if (j == 2) {
			k = 1;
			l = 0;
			m = 3;
		}
		if (j == 3) {
			k = 1;
			l = 2;
			m = 0;
		}
		area = this->fastSphericalTriangleArea(
		        cell->vertex(j)->point(), cell->vertex(k)->point().point(), cell->vertex(l)->point().point(), cell->vertex(m)->point().point());
		return area;
	}

	template <class Tesselation> std::vector<Real> FlowBoundingSphere<Tesselation>::getCellVelocity(Real X, Real Y, Real Z)
	{
		RTriangulation&   Tri = T[noCache ? (!currentTes) : currentTes].Triangulation();
		CellHandle        cell = Tri.locate(CGT::Sphere(X, Y, Z));
		std::vector<Real> velocityVector { cell->info().averageVelocity()[0], cell->info().averageVelocity()[1], cell->info().averageVelocity()[2] };
		return velocityVector;
	}

	template <class Tesselation> Real FlowBoundingSphere<Tesselation>::getCellVolume(Real X, Real Y, Real Z)
	{
		if (noCache && T[!currentTes].Max_id() <= 0) return 0; //the engine never solved anything
		RTriangulation& Tri = T[noCache ? (!currentTes) : currentTes].Triangulation();
		CellHandle      cell = Tri.locate(CGT::Sphere(X, Y, Z));
		return cell->info().volume(); //sqrt( cell->info().averageVelocity().squared_length());
	}


} //namespace CGT

} // namespace yade

#endif //FLOW_ENGINE
