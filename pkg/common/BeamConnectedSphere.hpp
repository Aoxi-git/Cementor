/*************************************************************************
*  Carlos Andrés del Valle Urberuaga
*  cdelv@unal.edu.co
*  2023
*************************************************************************/
/*
Implementation of BeamConnectedSphere class. This class is a child of the Sphere class, 
and it is used to connect two spheres with a beam. The beam uses a FEM-like formulation
to calculate the forces between the two spheres. The interaction has 6 degrees of freedom,
therefore the deformation and bending of the beam are considered. 
The beam class stores the necessary information for interaction.
The beams are elastic but support Rayleigh damping. 
The beams can also break according to the Mohr-Coulomb failure criterion.
The law2 that handles the interaction is also implemented. 
When two BeamConnectedSpheres are in contact but, no beam is present,
Law2_ScGeom_MindlinPhys_Mindlin is used. By default, the beam is created assumed to have a
circular cross-section. However, the user can change the beam geometry 
parameters to support different geometries.
*/
#pragma once
#include <core/Body.hpp>
#include <pkg/dem/ScGeom.hpp>
#include <pkg/common/Sphere.hpp>
#include <pkg/dem/HertzMindlin.hpp>

namespace yade { // Cannot have #include directive inside.

//!##################	BeamConnectedSphere  #####################
class BeamConnectedSphere : public Sphere {
public:
	virtual ~BeamConnectedSphere(){};
	void  addConnection(shared_ptr<Body> body);
	void  delConnection(int id);
    vector<shared_ptr<Body>> getConnections() const {return ConnList;};
	bool isConnectedTo(shared_ptr<Body> body);

	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS_CTOR_PY(BeamConnectedSphere, Sphere,"The class represents a sphere that can be connected to other beam-connected spheres with beams. The connections are stored in the ConnList attribute.",
		((vector<shared_ptr<Body>>,ConnList,,Attr::hidden,"List of :yref:`Beams` the sphere is connected to.")),
		/*ctor*/
		createIndex();,
		/*py*/
		.def("addConnection",&BeamConnectedSphere::addConnection,(boost::python::arg("Body")),"Adds a Beam to the connection list.")
		.def("delConnection",&BeamConnectedSphere::delConnection,(boost::python::arg("id")),"Removes a Beam from the connection list. It also removes the connection from the other node and eliminates the reference to the nodes that the beam stores.")
		.def("getConnections",&BeamConnectedSphere::getConnections,"Returns the list of connected :yref:`Beam`s.") 
		.def("isConnectedTo",&BeamConnectedSphere::isConnectedTo,(boost::python::arg("Body")),"Returns True if the sphere is connected to the Body passed as argument.")
	);
	// clang-format on
	REGISTER_CLASS_INDEX(BeamConnectedSphere, Sphere);
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(BeamConnectedSphere);

//!##################	BEAM   #####################
class Beam : public Shape {
public:
	virtual ~Beam(){};            
	void setNodes(shared_ptr<Body> Node1, shared_ptr<Body> Node2);
	void setDimensions(Real Radius, Real Length);
	void setMaterialProperties(Real Young_modulus, Real Shear_modulus, Real Density_);
	void setRayleighDampingCoefficients(Real A0, Real A1);
	void setBeamGeometry(Real A_, Real L_, Real Ix_, Real Iy_, Real Iz_, Real J_);
	void setSectionModulus(Real W);
	vector<Real> getEingenvaluesList(void);
	MatrixXr getStiffnessMatrix(void); 
	MatrixXr getMassMatrix(void);
	MatrixXr getDampingMatrix(void);
	void configureBeam(shared_ptr<Body> Node1, shared_ptr<Body> Node2);
	bool isFractured(Vector3r Force, Vector3r Torque);

	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS_CTOR_PY(Beam,Shape,"Beam class. Stores all the necessary information for interaction between two :yref:`BeamConnectedSphere`s. The beam is always aligned with the x-axis (in the beam reference frame) and has a circular cross-section by default. The geometric properties can be changed to support different geometries.",
		((shared_ptr<Body>, node1, nullptr,Attr::readonly,"First :yref:`Body` the Beam is connected to.")) 
		((shared_ptr<Body>, node2, nullptr,Attr::readonly,"Second :yref:`Body` the Beam is connected to."))
		((Real, L0, 0.0,Attr::readonly,"Initial Length of the Beam."))
		((Real, radius, 0.0,Attr::readonly,"Beam cross-section radius."))
		((Real, A, 0.0,Attr::readonly,"Beam cross-section area."))
		((Real, Ix, 0.0,Attr::readonly,"Beam second moment of area in the x direction.")) 
		((Real, Iy, 0.0,Attr::readonly,"Beam second moment of area in y directions."))
		((Real, Iz, 0.0,Attr::readonly,"Beam second moment of area in z directions."))
		((Real, J, 0.0,Attr::readonly,"Beam torsional constant of the cross-section."))
		((Real, E, 0.0,Attr::readonly,"Beam Young modulus."))
		((Real, G, 0.0,Attr::readonly,"Beam Shear modulus."))
		((Real, Density, 0.0,Attr::readonly,"Beam density."))
		((Real, Damping, 0.0,,"Beam Rayleigh damping fraction."))
		((Real, a0, 0.0,Attr::readonly,"Rayleigh damping: Mass matrix coefficient."))
		((Real, a1, 0.0,Attr::readonly,"Rayleigh damping: Stiffness matrix coefficient."))
		((bool, fracture, false,,"Bool that says whether or not the Beam can break."))
		((Real, Phi, 0.0,,"The angle of internal friction for Mohr Culomb failure criteria."))
		((Real, Cohesion, 0.0,,"Material cohesion for Mohr Culomb failure criteria."))
		((Real, w, 0.0,Attr::readonly,"Beam section modulus."))
		((Quaternionr, BeamInitialOrientation, Quaternionr(1.0,0.0,0.0,0.0),Attr::readonly,"Beam initial orientation."))
		((Quaternionr, NodeInitialOrientation1, Quaternionr(1.0,0.0,0.0,0.0),Attr::readonly,"Node1 :yref:`Body` initial orientation.")) // no need to use geom in the law if we store the initial orientations
		((Quaternionr, NodeInitialOrientation2, Quaternionr(1.0,0.0,0.0,0.0),Attr::readonly,"Node2 :yref:`Body` initial orientation.")),
		createIndex();, /*ctor*/
		/*py*/
		.def("setNodes",&Beam::setNodes,(boost::python::arg("Body"),boost::python::arg("Body")),"Receives the two :yref:`Body` the beam is connected to, checks that they are BeamConnectedSpheres, stores a pointer to the bodies, the initial orientations, and calculates the beam's initial orientation.")
		.def("setBeamDimensions",&Beam::setDimensions,(boost::python::arg("Radius"),boost::python::arg("Length")),"Sets the beam dimensions and calculates the beam geometry dependant properties: A, Ix, I, J. The beam is assumed to be a cylinder. When using this function, check that the provided length is consistent with the distance between the nodes. Otherwise, abnormally big forces can be created.")
		.def("setMaterialProperties",&Beam::setMaterialProperties,(boost::python::arg("Young_modulus"),boost::python::arg("Shear_modulus"),boost::python::arg("Density_")),"Sets the beam material properties and uses them to calculate the corresponding Rayleigh damping coefficients.")
		.def("setRayleighDampingCoefficients",&Beam::setRayleighDampingCoefficients,(boost::python::arg("A0"),boost::python::arg("A1")),"Sets the Rayleigh damping coefficients.")
		.def("setBeamGeometry",&Beam::setBeamGeometry,(boost::python::arg("A_"),boost::python::arg("L_"),boost::python::arg("Ix_"),boost::python::arg("Iy_"),boost::python::arg("Iz_"),boost::python::arg("J_")),"Sets the beam geometry properties. This way the user can create a custom beam. When using this function, check that the provided length is consistent with the distance between the nodes. Otherwise, abnormally big forces can be created.")
		.def("setSectionModulus",&Beam::setSectionModulus,(boost::python::arg("W")),"Sets the beam section modulus.")
		.def("getEingenvaluesList",&Beam::getEingenvaluesList,"Returns a list of the beam oscillation eigenmodes. They are the analytical solution to the eigenvalue problem det(K - omega^2 M) = 0. The function returns omega. There are only 6 eigenmodes because the other 6 are 0.")
		.def("getStiffnessMatrix",&Beam::getStiffnessMatrix,"Returns the beam stiffness matrix (K).")
		.def("getMassMatrix",&Beam::getMassMatrix,"Returns the beam mass matrix (M).") 
		.def("getDampingMatrix",&Beam::getDampingMatrix,"Returns the beam Rayleigh damping matrix (D = a0 M + a1 K).")
		.def("configureBeam",&Beam::configureBeam,(boost::python::arg("Body"),boost::python::arg("Body")),"Configures the beam using the properties of the nodes. The beam radius is the minimum of the radiuses of the nodes, the length is the distance between the nodes, and the material properties are the harmonic average between the properties of the nodes.") 
		.def("isFractured",&Beam::isFractured,(boost::python::arg("Force"),boost::python::arg("Torque")),"Returns true if the beam will fracture under the given forces and torques acting on it. The failure criteria is Mohr-Coulomb. The force and torque have to be in the beam reference frame.")
	);
	// clang-format on
	REGISTER_CLASS_INDEX(Beam, Shape);
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(Beam);

class Law2_ScGeom_MindlinPhys_BeamConnectedSphere : public LawFunctor {
public:
	bool go(shared_ptr<IGeom>& _geom, shared_ptr<IPhys>& _phys, Interaction* I) override;
	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS_CTOR_PY(Law2_ScGeom_MindlinPhys_BeamConnectedSphere, LawFunctor,"Law between a frictional :yref:`BeamConnectedSphere` and a frictional :yref:`BeamConnectedSphere`. When there's no beam connection, Law2_ScGeom_MindlinPhys_Mindlin is called.",
		((bool,neverErase,true,,"Keep interactions even if particles go away from each other."))
		((bool,includeAdhesion,false,,"bool to include the adhesion force following the DMT formulation. If true, also the normal elastic energy takes into account the adhesion effect."))
		((bool,calcEnergy,false,,"bool to calculate energy terms (shear potential energy, dissipation of energy due to friction and dissipation of energy due to normal and tangential damping)"))
		((bool,includeMoment,false,,"bool to consider rolling resistance (if :yref:`Ip2_FrictMat_FrictMat_MindlinPhys::eta` is 0.0, no plastic condition is applied.)"))
		((shared_ptr<Law2_ScGeom_MindlinPhys_Mindlin>, Hertz, nullptr, Attr::hidden, "Hertz law for the interaction between non connected bodies.")),
		Hertz.reset(new Law2_ScGeom_MindlinPhys_Mindlin);
		Hertz->includeAdhesion = includeAdhesion;
		Hertz->calcEnergy = calcEnergy;
		Hertz->includeMoment = includeMoment;
		timingDeltas=shared_ptr<TimingDeltas>(new TimingDeltas);,
	);
	// clang-format on
	FUNCTOR2D(ScGeom, MindlinPhys); // Same as Law2_ScGeom_MindlinPhys_Mindlin for compatibility
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(Law2_ScGeom_MindlinPhys_BeamConnectedSphere);

} // namespace yade