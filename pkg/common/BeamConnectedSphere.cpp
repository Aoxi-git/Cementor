/*************************************************************************
*  Carlos Andrés del Valle Urberuaga
*  cdelv@unal.edu.co
*  2023
*************************************************************************/
#include <pkg/common/BeamConnectedSphere.hpp>
#include <algorithm>
#include <iostream>
#include <lib/base/Math.hpp>
#include <lib/high-precision/Constants.hpp>

namespace yade { // Cannot have #include directive inside.

//!##################	BeamConnectedSphere  #####################
YADE_PLUGIN((BeamConnectedSphere));
CREATE_LOGGER(BeamConnectedSphere);

void BeamConnectedSphere::addConnection(shared_ptr<Body> body) {
    // Check that the added object is a BeamConnectedSphere
    if (body->shape->getClassIndex() != Beam::getClassIndexStatic()){
        LOG_FATAL("Trying to add a non-Beam object to the list of connections.")
        throw std::runtime_error("BeamConnectedSphere::addConnection: Trying to add a non-Beam object to the list of connections.");
    }

    // Check that the added object is not itself
    if (body->shape == shared_from_this()){
        LOG_FATAL("Trying to add itself to the list of connections.")
        throw std::runtime_error("BeamConnectedSphere::addConnection: Trying to add itself to the list of connections.");
    }

    // Check that the added object is not already in the list
    if (std::find(ConnList.begin(), ConnList.end(), body) != ConnList.end())
        LOG_WARN("Adding an already existing connection.");
        
    // Add the new beam to the list
    ConnList.push_back(body);

    // Sort the list
    std::sort(ConnList.begin(), ConnList.end());

    // Make sure there are no duplicates
    ConnList.erase(std::unique(ConnList.begin(), ConnList.end()), ConnList.end());
}

void BeamConnectedSphere::delConnection(int id) {
    // Check that the id is valid
    if (id < 0 || static_cast<unsigned int>(id) >= ConnList.size()){
        LOG_FATAL("Trying to delete a non-existing connection.");
        throw std::runtime_error("BeamConnectedSphere::delConnection: Trying to delete a non-existing connection.");
    }
    // Get the beam from the ConnList to be deleted
    shared_ptr<Beam> del_beam = YADE_PTR_CAST<Beam>(ConnList[id]->shape);

    // Get the nodes of the beam connection
    shared_ptr<BeamConnectedSphere> node1 = YADE_PTR_CAST<BeamConnectedSphere>(del_beam->node1->shape);
    shared_ptr<BeamConnectedSphere> node2 = YADE_PTR_CAST<BeamConnectedSphere>(del_beam->node2->shape);

    // Remove the beam from the connection list of the two nodes
    node1->ConnList.erase(std::remove(node1->ConnList.begin(), node1->ConnList.end(), ConnList[id]), node1->ConnList.end());
    node2->ConnList.erase(std::remove(node2->ConnList.begin(), node2->ConnList.end(), ConnList[id]), node2->ConnList.end());

    // Remove the pointer to the nodes of the beam to be deleted
    del_beam->node1 = nullptr;
    del_beam->node2 = nullptr;
}

bool BeamConnectedSphere::isConnectedTo(shared_ptr<Body> body){
    // Check that the body is a BeamConnectedSphere
    if (body->shape->getClassIndex() != getClassIndex()){
        LOG_WARN("Trying to check if a non-BeamConnectedSphere object is connected to this node.")
        return false;
    }

    // Check that the body is not it self
    if (body->shape == shared_from_this()){
        LOG_WARN("Trying to check if a body is connected to itself.")
        return false;
    }

    // Find common elements in the two lists.
    const vector<shared_ptr<Body>> &bodyConnections = YADE_PTR_CAST<BeamConnectedSphere>(body->shape)->getConnections();
    for (const auto &i : bodyConnections)
        if (std::find(ConnList.begin(), ConnList.end(), i) != ConnList.end())
            return true;

    return false;
}

//!##################	BEAM   #####################
YADE_PLUGIN((Beam));
CREATE_LOGGER(Beam);

void Beam::setNodes(shared_ptr<Body> Node1, shared_ptr<Body> Node2){
    // Check that the nodes are defined
    if (!Node1 || !Node2){
        LOG_FATAL("Node1 or Node2 is not defined. They have to be set before calling this function.");
        throw std::runtime_error("Beam::setNodes: Node1 or Node2 is not defined. They have to be set before calling this function.");
    }

    // Check that the nodes are BeamConnectedSphere
    if (Node1->shape->getClassIndex() != BeamConnectedSphere::getClassIndexStatic() || Node2->shape->getClassIndex() != BeamConnectedSphere::getClassIndexStatic()){
        LOG_FATAL("Node1 or Node2 is not a BeamConnectedSphere.");
        throw std::runtime_error("Beam::setNodes: Node1 or Node2 is not a BeamConnectedSphere.");
    }

    // Check that the nodes are not the same
    if (Node1 == Node2){
        LOG_FATAL("Node1 and Node2 are the same.");
        throw std::runtime_error("Beam::setNodes: Node1 and Node2 are the same.");
    }

    // Store the pointers to the nodes
    node1 = Node1;
    node2 = Node2;

    // Calculate and store beam initial orientation
    BeamInitialOrientation = Quaternionr::FromTwoVectors(Vector3r(1.0, 0.0, 0.0), (node2->state->pos - node1->state->pos).normalized());
    BeamInitialOrientation.normalize();
    
    // Store beam nodes initial orientation
    NodeInitialOrientation1 = node1->state->ori.normalized();
    NodeInitialOrientation2 = node2->state->ori.normalized();
}

void Beam::setDimensions(Real Radius, Real Length){
    if (Radius <= 0.0 || Length <= 0.0) {
        LOG_FATAL("Radius or Length values are not valid. They have to be positive and != 0 numbers");
        throw std::runtime_error("Beam::setDimensions: Radius or Length values are not valid. They have to be positive and != 0 numbers");
    }

    // Store beam dimensions
    radius = Radius;
    L0 = Length;

    // Calculate beam geometric properties
    A = Mathr::PI*radius*radius;                     // Cross-sectional area
    Ix = 2.0*L0*A*sqrt(A/Mathr::PI)/(3.0*Mathr::PI); // Second moment of area in dir x
    Iy = A*A/(4.0*Mathr::PI);                         // The second moment of area in dir y-z, as cross-section is assumed circular Iy = Iz = I
    Iz = A*A/(4.0*Mathr::PI);                         // The second moment of area in dir y-z, as cross-section is assumed circular Iy = Iz = I
    J = A*A/(2.0*Mathr::PI);                         // Beam torsional constant of the cross-section.

    // Calculate Section modulus (for failure criterion)
    w = Iy/radius; // I/r = Iy/y_max = Iz/z_max (circular cross-section)
}

void Beam::setMaterialProperties(Real Young_modulus, Real Shear_modulus, Real Density_){
    if (Young_modulus <= 0.0 || Shear_modulus <= 0.0 || Density_ <= 0.0) {
        LOG_FATAL("Young_modulus, Shear_modulus, or Density values are not valid. They have to be positive and != 0 numbers.");
        throw std::runtime_error("Beam::setMaterialProperties: Young_modulus, Shear_modulus, or Density values are not valid. They have to be positive and != 0 numbers.");
    }

    // Check Poisson ratio
    Real aux_poisson = Young_modulus/(2.0*Shear_modulus) - 1.0;
    if (aux_poisson < -1.0 || aux_poisson > 0.5)
        LOG_WARN("Poisson ratio should be between -1 and 0.5. The properties provided give the value: " << aux_poisson << " .");

    // Store beam material properties
    E = Young_modulus;
    G = Shear_modulus;
    Density = Density_;

    // Rayleigh damping parameters
    // Solution to the eigenvalue problem det(K - ω²M) = 0.
    // As the matrices are 12x12 there should be 12 eigenvalues. However, 6 of them are 0.
    // Gets ω sorted from smallest to largest of the 6 non-zero eigenvalues.
    vector<Real> omega = getEingenvaluesList();

    // Select normal modes for the damping matrix. I use the first 2 modes.
    Real w1 = omega[0];
    Real w2 = omega[1];

    // Compute Rayleigh damping coefficients
    a0 = 2.0*Damping*w1*w2/(w1+w2);
    a1 = 2.0*Damping/(w1+w2);
}

void Beam::setRayleighDampingCoefficients(Real A0, Real A1){ 
    if (A0 < 0.0 || A1 < 0.0) 
        LOG_WARN("A0 or A1 value is negative.");
 
    // Store Rayleigh damping coefficients
    a0 = A0;
    a1 = A1;
}

void Beam::setBeamGeometry(Real A_, Real L_, Real Ix_, Real Iy_, Real Iz_, Real J_){
    if (A_ <= 0.0 || L_ <= 0.0 || Ix_ <= 0.0 || Iy_ <= 0.0 || Iz_ <= 0.0 || J_ <= 0.0) {
        LOG_FATAL("A, L, Ix, Iy, Iz, or J values are not valid. They have to be positive and != 0 numbers.");
        throw std::runtime_error("Beam::setBeamGeometry: A, L, Ix, Iy, Iz, or J values are not valid. They have to be positive and != 0 numbers.");
    }

    // Store beam dimensions
    A = A_;
    L0 = L_;
    Ix = Ix_;
    Iy = Iy_;
    Iz = Iz_;
    J = J_;
}

void Beam::setSectionModulus(Real W){
    if (W <= 0.0) {
        LOG_FATAL("W value is not valid. It has to be a positive and != 0 number.");
        throw std::runtime_error("Beam::setSectionModulus: W value is not valid. It has to be a positive and != 0 number.");
    }

    // Store section modulus
    w = W;
}

vector<Real> Beam::getEingenvaluesList(void){
    // Analytical solution to the eigenvalue problem det(K - ω²M) = 0.
    // As the matrices are 12x12, there should be 12 eigenvalues. However, 6 of them are 0.
    vector<Real> omega2 = { 
        720.0*E*Iz/(A*pow(L0,4)*Density), 
        8400.0*E*Iz/(A*pow(L0,4)*Density), 
        720.0*E*Iy/(A*pow(L0,4)*Density), 
        8400.0*E*Iy/(A*pow(L0,4)*Density), 
        12.0*G*J/(Ix*pow(L0,2)*Density), 
        12.0*E/(pow(L0,2)*Density)
    };
    // Calculate the square root of the eigenvalues
    std::transform(omega2.begin(), omega2.end(), omega2.begin(), [](Real val){return sqrt(val);});

    // Sort the eigenvalues from smallest to largest
    std::sort(omega2.begin(), omega2.end());
    
    return omega2;
}
void Beam::configureBeam(shared_ptr<Body> Node1, shared_ptr<Body> Node2){
    // Add beam nodes. 
    setNodes(Node1, Node2);

    // Add beam geometry properties. We assume that the beam cross-section is a circle
    setDimensions(math::min(YADE_PTR_CAST<BeamConnectedSphere>(node1->shape)->radius, YADE_PTR_CAST<BeamConnectedSphere>(node2->shape)->radius), (node2->state->pos - node1->state->pos).norm());

    // Add beam material properties
    Real E1 = node1->material->cast<FrictMat>().young; // Use cast to do runtime type checking
    Real E2 = node1->material->cast<FrictMat>().young;
    Real nu1 = node1->material->cast<FrictMat>().poisson;
    Real nu2 = node2->material->cast<FrictMat>().poisson;
    Real G1 = E1/(2.0*(1.0+nu1));
    Real G2 = E2/(2.0*(1.0+nu2));
    Real rho1 = node1->material->density;
    Real rho2 = node2->material->density;
    setMaterialProperties(2.0*E1*E2/(E1+E2), 2.0*G1*G2/(G1+G2), 2.0*rho1*rho2/(rho1+rho2));
}

MatrixXr Beam::getStiffnessMatrix(void){
    /* Matrix{12, 12}(
        L*L*A ,      0.0,       0.0,        0.0,        0.0,        0.0, -L*L*A,       0.0,      0.0,        0.0,        0.0,        0.0,
        0.0   ,  12.0*Iz,       0.0,        0.0,        0.0,   6.0*L*Iz,    0.0,  -12.0*Iz,      0.0,        0.0,        0.0,   6.0*L*Iz,
        0.0   ,      0.0,   12.0*Iy,        0.0,  -6.0*L*Iy,        0.0,    0.0,       0.0, -12.0*Iy,        0.0,  -6.0*L*Iy,        0.0,
        0.0   ,      0.0,       0.0,  G*L*L*J/E,        0.0,        0.0,    0.0,       0.0,      0.0, -G*L*L*J/E,        0.0,        0.0,
        0.0   ,      0.0, -6.0*L*Iy,        0.0, 4.0*L*L*Iy,        0.0,    0.0,       0.0, 6.0*L*Iy,        0.0, 2.0*L*L*Iy,        0.0,
        0.0   , 6.0*L*Iz,       0.0,        0.0,        0.0, 4.0*L*L*Iz,    0.0, -6.0*L*Iz,      0.0,        0.0,        0.0, 2.0*L*L*Iz,
        -L*L*A,      0.0,       0.0,        0.0,        0.0,        0.0,  L*L*A,       0.0,      0.0,        0.0,        0.0,        0.0,
        0.0   , -12.0*Iz,       0.0,        0.0,        0.0,  -6.0*L*Iz,    0.0,   12.0*Iz,      0.0,        0.0,        0.0,  -6.0*L*Iz,
        0.0   ,      0.0,  -12.0*Iy,        0.0,   6.0*L*Iy,        0.0,    0.0,       0.0,  12.0*Iy,        0.0,   6.0*L*Iy,        0.0,
        0.0   ,      0.0,       0.0, -G*L*L*J/E,        0.0,        0.0,    0.0,       0.0,      0.0,  G*L*L*J/E,        0.0,        0.0,
        0.0   ,      0.0, -6.0*L*Iy,        0.0, 2.0*L*L*Iy,        0.0,    0.0,       0.0, 6.0*L*Iy,        0.0, 4.0*L*L*Iy,        0.0,
        0.0   , 6.0*L*Iz,       0.0,        0.0,        0.0, 2.0*L*L*Iz,    0.0, -6.0*L*Iz,      0.0,        0.0,        0.0, 4.0*L*L*Iz
    )*(E/L^3) */
    MatrixXr K = MatrixXr::Zero(12, 12);
    K(0, 0) = L0*L0*A;
    K(0, 6) = -L0*L0*A;
    K(1, 1) = 12.0*Iz;
    K(1, 5) = 6.0*L0*Iz;
    K(1, 7) = -12.0*Iz;
    K(1, 11) = 6.0*L0*Iz;
    K(2, 2) = 12.0*Iy;
    K(2, 4) = -6.0*L0*Iy;
    K(2, 8) = -12.0*Iy;
    K(2, 10) = -6.0*L0*Iy;
    K(3, 3) = G*L0*L0*J/E;
    K(3, 9) = -G*L0*L0*J/E;
    K(4, 2) = -6.0*L0*Iy;
    K(4, 4) = 4.0*L0*L0*Iy;
    K(4, 8) = 6.0*L0*Iy; 
    K(4, 10) = 2.0*L0*L0*Iy; 
    K(5, 1) = 6.0*L0*Iz;
    K(5, 5) = 4.0*L0*L0*Iz;
    K(5, 7) = -6.0*L0*Iz; 
    K(5, 11) = 2.0*L0*L0*Iz; 
    K(6, 0) = -L0*L0*A;
    K(6, 6) = L0*L0*A;
    K(7, 1) = -12.0*Iz;
    K(7, 5) = -6.0*L0*Iz;
    K(7, 7) = 12.0*Iz;
    K(7, 11) = -6.0*L0*Iz;
    K(8, 2) = -12.0*Iy;
    K(8, 4) = 6.0*L0*Iy;
    K(8, 8) = 12.0*Iy;
    K(8, 10) = 6.0*L0*Iy;
    K(9, 9) = G*L0*L0*J/E;
    K(9, 3) = -G*L0*L0*J/E; 
    K(10, 2) = -6.0*L0*Iy;
    K(10, 4) = 2.0*L0*L0*Iy; 
    K(10, 8) = 6.0*L0*Iy;
    K(10, 10) = 4.0*L0*L0*Iy;
    K(11, 1) = 6.0*L0*Iz;
    K(11, 5) = 2.0*L0*L0*Iz;
    K(11, 7) = -6.0*L0*Iz;
    K(11, 11) = 4.0*L0*L0*Iz;
    K *= E/pow(L0, 3);
    return K;
}

MatrixXr Beam::getMassMatrix(void){
    /*Matrix{12, 12}(
        140.0,     0.0,     0.0,      0.0,      0.0,      0.0,  70.0,     0.0,     0.0,      0.0,      0.0,      0.0,
        0.0  ,   156.0,     0.0,      0.0,      0.0,   22.0*L,   0.0,    54.0,     0.0,      0.0,      0.0,  -13.0*L,
        0.0  ,     0.0,   156.0,      0.0,  -22.0*L,      0.0,   0.0,     0.0,    54.0,      0.0,   13.0*L,      0.0,
        0.0  ,     0.0,     0.0, 140.0*r2,      0.0,      0.0,   0.0,     0.0,     0.0,  70.0*r2,      0.0,      0.0,
        0.0  ,     0.0, -22.0*L,      0.0,  4.0*L*L,      0.0,   0.0,     0.0, -13.0*L,      0.0, -3.0*L*L,      0.0,
        0.0  ,  22.0*L,     0.0,      0.0,      0.0,  4.0*L*L,   0.0,  13.0*L,     0.0,      0.0,      0.0,  -3.0*L*L,
        70.0 ,     0.0,     0.0,      0.0,      0.0,      0.0, 140.0,     0.0,     0.0,      0.0,      0.0,      0.0,
        0.0  ,    54.0,     0.0,      0.0,      0.0,   13.0*L,   0.0,   156.0,     0.0,      0.0,      0.0,  -22.0*L,
        0.0  ,     0.0,    54.0,      0.0,  -13.0*L,      0.0,   0.0,     0.0,   156.0,      0.0,   22.0*L,      0.0,
        0.0  ,     0.0,     0.0,  70.0*r2,      0.0,      0.0,   0.0,     0.0,     0.0, 140.0*r2,      0.0,      0.0,
        0.0  ,     0.0,  13.0*L,      0.0, -3.0*L*L,      0.0,   0.0,     0.0,  22.0*L,      0.0,  4.0*L*L,      0.0,
        0.0  , -13.0*L,     0.0,      0.0,      0.0, -3.0*L*L,   0.0, -22.0*L,     0.0,      0.0,      0.0,  4.0*L*L,
    )*(ρ*A*L/420.0)*/
    Real r2 = Ix/A; // Radius of gyration of the cross-section 

    MatrixXr M = MatrixXr::Zero(12, 12);
    M(0, 0) = 140.0;
    M(0, 6) = 70.0;
    M(1, 1) = 156.0;
    M(1, 5) = 22.0*L0;
    M(1, 7) = 54.0;
    M(1, 11) = -13.0*L0;
    M(2, 2) = 156.0;
    M(2, 4) = -22.0*L0;
    M(2, 8) = 54.0;
    M(2, 10) = 13.0*L0;
    M(3, 3) = 140.0*r2;
    M(3, 9) = 70.0*r2;
    M(4, 2) = -22.0*L0;
    M(4, 4) = 4.0*L0*L0;
    M(4, 8) = -13.0*L0;
    M(4, 10) = -3.0*L0*L0;
    M(5, 1) = 22.0*L0;
    M(5, 5) = 4.0*L0*L0;
    M(5, 7) = 13.0*L0;
    M(5, 11) = -3.0*L0*L0;
    M(6, 0) = 70.0;
    M(6, 6) = 140.0;
    M(7, 1) = 54.0;
    M(7, 5) = 13.0*L0;
    M(7, 7) = 156.0;
    M(7, 11) = -22.0*L0;
    M(8, 2) = 54.0;
    M(8, 4) = -13.0*L0;
    M(8, 8) = 156.0;
    M(8, 10) = 22.0*L0;
    M(9, 3) = 70.0*r2;
    M(9, 9) = 140.0*r2;
    M(10, 2) = 13.0*L0;
    M(10, 4) = -3.0*L0*L0;
    M(10, 8) = 22.0*L0;
    M(10, 10) = 4.0*L0*L0;
    M(11, 1) = -13.0*L0;
    M(11, 5) = -3.0*L0*L0;
    M(11, 7) = -22.0*L0;
    M(11, 11) = 4.0*L0*L0;
    M *= (Density*A*L0/420.0);
    return M;
}

MatrixXr Beam::getDampingMatrix(void){
    return a0*getMassMatrix() + a1*getStiffnessMatrix(); // Rayleigh damping
}

bool Beam::isFractured(Vector3r Force, Vector3r Torque){
    // Calculate stresses in the beam
    Real sigma_x = Force(0)/A + Torque(1)/w + Torque(2)/w; // Fx/A + My/w + Mz/w
    Real tau_y = 4.0*Force(1)/(3.0*A);                     // 4*Fy/(3*A)
    Real tau_z = 4.0*Force(2)/(3.0*A) + Torque(0)/w;       // 4*Fz/(3*A) + Mx/w  // The 4/3 factor is because the beam has a circular cross-section

    // Calculate the principal stress
    Real sigma_2 = sqrt(0.25*sigma_x*sigma_x + tau_y*tau_y + tau_z*tau_z);
    Real sigma_1 = 0.5*sigma_x + sigma_2;
    Real sigma_3 = 0.5*sigma_x - sigma_2;

    // Mohr Culomb stresses
    Real sigma = 0.5*(sigma_1 + sigma_3) + 0.5*(sigma_1 - sigma_3)*sin(Phi);
    Real tau = 0.5*(sigma_1 - sigma_3)*cos(Phi);

    // Check if the beam is fractured
    // the sing in σ*tan(ϕ) means that compression is negative. The Tensile is positive. 
    // This is a convention and can be flipped by changing the sign.
    if (math::abs(tau) >= Cohesion - sigma*tan(Phi))
        return true;

    return false;
}

//!##################	LAW   #####################
YADE_PLUGIN((Law2_ScGeom_MindlinPhys_BeamConnectedSphere));
CREATE_LOGGER(Law2_ScGeom_MindlinPhys_BeamConnectedSphere);

bool Law2_ScGeom_MindlinPhys_BeamConnectedSphere::go(shared_ptr<IGeom>& ig, shared_ptr<IPhys>& ip, Interaction* contact)
{
    // Get the BeamConnectedSpheres ids
    const int id1 = contact->getId1(), id2 = contact->getId2();
    Hertz->scene = scene; // Set the scene for the Hertz law. Cant do it in the constructor because the scene is not yet created. This is faster than checking if the scene is set every time the law is called.

    // Check if the bodies are BeamConnectedSpheres, if not, call Hertz law
    if (Body::byId(id1)->shape->getClassIndex() != BeamConnectedSphere::getClassIndexStatic() || Body::byId(id2)->shape->getClassIndex() != BeamConnectedSphere::getClassIndexStatic())
        return Hertz->go(ig, ip, contact);
    
    // Check if they are connected to the same beam
    int Beam_id = -1; // Id of connecting beam in the connection list of the BeamConnectedSphere(id2)
    shared_ptr<Beam> CurrentBeam; { // Create a new scope to delete the lists after use
        vector<shared_ptr<Body>> list2 = YADE_PTR_CAST<BeamConnectedSphere>((Body::byId(id2)->shape))->getConnections();
        for (const auto &i : YADE_PTR_CAST<BeamConnectedSphere>((Body::byId(id1)->shape))->getConnections()){ // Find common elements in the two lists.
            const auto it = std::find(list2.begin(), list2.end(), i);
            if ( it != list2.end()){
                CurrentBeam = YADE_PTR_CAST<Beam>(i->shape);
                Beam_id = it - list2.begin();
                break; // Only one beam can connect the two bodies
            }
        }
    } // Close scope

    // If no beam is connected, use Hertz law.
    if (Beam_id < 0)
        return Hertz->go(ig, ip, contact);

    // Calculate the orientation relative to the beam of the two nodes
    // This is how much each node has rotated + the initial orientation of the beam
    // The difference between the two orientations is the twisting of the beam
    // This is different from what SCgeom6D does because it does not account for the initial orientation of the beam
    Quaternionr Dq1 = (Body::byId(id1, scene)->state->ori*CurrentBeam->NodeInitialOrientation1.conjugate()*CurrentBeam->BeamInitialOrientation);
    Quaternionr Dq2 = (Body::byId(id2, scene)->state->ori*CurrentBeam->NodeInitialOrientation2.conjugate()*CurrentBeam->BeamInitialOrientation);
    Quaternionr Beam_q = Dq1.slerp(0.5, Dq2);       // Set the beam orientation to the average of the two (they should be very similar, although not always)
    Beam_q.normalize();                            // Normalize the quaternion. Eigen does not do it inside SLERP
    const Matrix3r R = Beam_q.toRotationMatrix(); // Beam - global transformation. Matrix rotation is faster than quaternion rotation

    // Account for periodic boundary conditions
    const Vector3r shift    = scene->isPeriodic ? scene->cell->intrShiftPos(contact->cellDist) : Vector3r::Zero();
	const Vector3r shiftVel = scene->isPeriodic ? scene->cell->intrShiftVel(contact->cellDist) : Vector3r::Zero();

    // Calculate the displacement vector in the beam frame
    // remember that the beam is always aligned with the X-axis in the beam coordinate system
    // We need to move to the frame at the center of the beam. That's why the division by 2
    const Vector3r Dx = 0.5*(R.transpose()*(Body::byId(id2, scene)->state->pos - Body::byId(id1, scene)->state->pos + shift) - Vector3r(CurrentBeam->L0,0.0,0.0));

    // Calculate the beam twisting and bending
    // We need to move to the frame at the center of the beam. That's why the division by 2
    AngleAxisr bending((Dq1.conjugate()*Dq2).normalized());
    if (bending.angle() > Mathr::PI) bending.angle() -= Mathr::TWO_PI; // Angle is between 0 and 2*pi, but should be between -pi and pi
    const Vector3r Dq = 0.5*bending.angle()*bending.axis();

    // Calculate relative velocity and angular velocity (for the Rayleigh damping)
    const Vector3r RelVel = 0.5*(R.transpose()*(Body::byId(id2, scene)->state->vel - Body::byId(id1, scene)->state->vel + shiftVel)); // Relative velocity in the beam coordinate system
    const Vector3r RelAngVel = 0.5*(R.transpose()*(Body::byId(id2, scene)->state->angVel - Body::byId(id1, scene)->state->angVel)); // Relative angular velocity in the beam coordinate system
    
    // Declare the force and torque vectors
    Vector3r force, torque1, torque2;
    /*
    // Force Calculated using Matrices. I leave this here because it could be useful in the future for incorporating plasticity
    // Also, it is much easier to understand the force model by looking at this than the obfuscated analytical solution I use below.
    // I use the analytical solution because it removes the need to store the matrices. It's faster and more accurate.
    // Commenting the analytical force, and uncommenting this will give the same results.   MatrixXr X = MatrixXr::Zero(12, 1);
    // Create the 12D dof vectors for the force calculation
    MatrixXr U = MatrixXr::Zero(12, 1);
    // Fill the dof vector (position and twisting)
    X(0,0)  = Dx(0);
    X(1,0)  = Dx(1);
    X(2,0)  = Dx(2);
    X(3,0)  = Dq(0);
    X(4,0)  = Dq(1);
    X(5,0)  = Dq(2);
    X(6,0)  = -Dx(0);
    X(7,0)  = -Dx(1);
    X(8,0)  = -Dx(2);
    X(9,0)  = -Dq(0);
    X(10,0) = -Dq(1);
    X(11,0) = -Dq(2);
    // Fill the dof vector (velocity and angular velocity)
    U(0,0)  = RelVel(0);
    U(1,0)  = RelVel(1);
    U(2,0)  = RelVel(2);
    U(3,0)  = RelAngVel(0);
    U(4,0)  = RelAngVel(1);
    U(5,0)  = RelAngVel(2);
    U(6,0)  = -RelVel(0);
    U(7,0)  = -RelVel(1);
    U(8,0)  = -RelVel(2);
    U(9,0)  = -RelAngVel(0);
    U(10,0) = -RelAngVel(1);
    U(11,0) = -RelAngVel(2);
    // Calculate the force 12D dof vector
    MatrixXr F = CurrentBeam->getStiffnessMatrix()*X + CurrentBeam->getDampingMatrix()*U;
    // Force
    force = Vector3r(F(0,0), F(1,0), F(2,0));
    // Torques are not just the negative of the other. Be careful with this!
    torque1 = Vector3r(F(3,0), F(4,0), F(5,0));
    torque2 = Vector3r(F(9,0), F(10,0), F(11,0));
    */
    // ANALYTICAL FORCE CALCULATION
    force = Vector3r( // Elastic force
        2.0*CurrentBeam->A*CurrentBeam->E*Dx(0)/CurrentBeam->L0,
        24.0*CurrentBeam->E*CurrentBeam->Iz*Dx(1)/pow(CurrentBeam->L0,3),
        24.0*CurrentBeam->E*CurrentBeam->Iy*Dx(2)/pow(CurrentBeam->L0,3)
    );
    torque1 = Vector3r( // Elastic torque
        2.0*CurrentBeam->G*CurrentBeam->J*Dq(0)/CurrentBeam->L0,
        2.0*CurrentBeam->E*CurrentBeam->Iy*Dq(1)/CurrentBeam->L0 - 12.0*CurrentBeam->E*CurrentBeam->Iy*Dx(2)/pow(CurrentBeam->L0,2),
        2.0*CurrentBeam->E*CurrentBeam->Iz*Dq(2)/CurrentBeam->L0 + 12.0*CurrentBeam->E*CurrentBeam->Iz*Dx(1)/pow(CurrentBeam->L0,2)
    );
    torque2 = Vector3r( // Elastic torque
        -2.0*CurrentBeam->G*CurrentBeam->J*Dq(0)/CurrentBeam->L0,
        -2.0*CurrentBeam->E*CurrentBeam->Iy*Dq(1)/CurrentBeam->L0 - 12.0*CurrentBeam->E*CurrentBeam->Iy*Dx(2)/pow(CurrentBeam->L0,2),
        -2.0*CurrentBeam->E*CurrentBeam->Iz*Dq(2)/CurrentBeam->L0 + 12.0*CurrentBeam->E*CurrentBeam->Iz*Dx(1)/pow(CurrentBeam->L0,2)
    );
    if (CurrentBeam->Damping > 0.0){
        force += Vector3r( // Damping force
            CurrentBeam->A*CurrentBeam->L0*CurrentBeam->a0*RelVel(0)*CurrentBeam->Density/6.0 + 2.0*CurrentBeam->A*CurrentBeam->E*CurrentBeam->a1*RelVel(0)/CurrentBeam->L0,
            CurrentBeam->A*pow(CurrentBeam->L0,2)*CurrentBeam->a0*RelAngVel(2)*CurrentBeam->Density/12.0 + 17.0*CurrentBeam->A*CurrentBeam->L0*CurrentBeam->a0*RelVel(1)*CurrentBeam->Density/70.0 + 24.0*CurrentBeam->E*CurrentBeam->Iz*CurrentBeam->a1*RelVel(1)/pow(CurrentBeam->L0,3),
            -CurrentBeam->A*pow(CurrentBeam->L0,2)*CurrentBeam->a0*RelAngVel(1)*CurrentBeam->Density/12.0 + 17.0*CurrentBeam->A*CurrentBeam->L0*CurrentBeam->a0*RelVel(2)*CurrentBeam->Density/70.0 + 24.0*CurrentBeam->E*CurrentBeam->Iy*CurrentBeam->a1*RelVel(2)/pow(CurrentBeam->L0,3)
        ); 
        torque1 += Vector3r( // Damping torque
            CurrentBeam->Ix*CurrentBeam->L0*CurrentBeam->a0*RelAngVel(0)*CurrentBeam->Density/6.0 + 2.0*CurrentBeam->G*CurrentBeam->J*CurrentBeam->a1*RelAngVel(0)/CurrentBeam->L0,
            CurrentBeam->A*pow(CurrentBeam->L0,3)*CurrentBeam->a0*RelAngVel(1)*CurrentBeam->Density/60.0 - 3.0*CurrentBeam->A*pow(CurrentBeam->L0,2)*CurrentBeam->a0*RelVel(2)*CurrentBeam->Density/140.0 + 2.0*CurrentBeam->E*CurrentBeam->Iy*CurrentBeam->a1*RelAngVel(1)/CurrentBeam->L0 - 12.0*CurrentBeam->E*CurrentBeam->Iy*CurrentBeam->a1*RelVel(2)/pow(CurrentBeam->L0,2),
            CurrentBeam->A*pow(CurrentBeam->L0,3)*CurrentBeam->a0*RelAngVel(2)*CurrentBeam->Density/60.0 + 3.0*CurrentBeam->A*pow(CurrentBeam->L0,2)*CurrentBeam->a0*RelVel(1)*CurrentBeam->Density/140.0 + 2.0*CurrentBeam->E*CurrentBeam->Iz*CurrentBeam->a1*RelAngVel(2)/CurrentBeam->L0 + 12.0*CurrentBeam->E*CurrentBeam->Iz*CurrentBeam->a1*RelVel(1)/pow(CurrentBeam->L0,2)
        );
        torque2 += Vector3r( // Damping torque
            -CurrentBeam->Ix*CurrentBeam->L0*CurrentBeam->a0*RelAngVel(0)*CurrentBeam->Density/6.0 - 2.0*CurrentBeam->G*CurrentBeam->J*CurrentBeam->a1*RelAngVel(0)/CurrentBeam->L0,
            -CurrentBeam->A*pow(CurrentBeam->L0,3)*CurrentBeam->a0*RelAngVel(1)*CurrentBeam->Density/60.0 - 3.0*CurrentBeam->A*pow(CurrentBeam->L0,2)*CurrentBeam->a0*RelVel(2)*CurrentBeam->Density/140 - 2.0*CurrentBeam->E*CurrentBeam->Iy*CurrentBeam->a1*RelAngVel(1)/CurrentBeam->L0 - 12.0*CurrentBeam->E*CurrentBeam->Iy*CurrentBeam->a1*RelVel(2)/pow(CurrentBeam->L0,2),
            -CurrentBeam->A*pow(CurrentBeam->L0,3)*CurrentBeam->a0*RelAngVel(2)*CurrentBeam->Density/60.0 + 3.0*CurrentBeam->A*pow(CurrentBeam->L0,2)*CurrentBeam->a0*RelVel(1)*CurrentBeam->Density/140 - 2.0*CurrentBeam->E*CurrentBeam->Iz*CurrentBeam->a1*RelAngVel(2)/CurrentBeam->L0 + 12.0*CurrentBeam->E*CurrentBeam->Iz*CurrentBeam->a1*RelVel(1)/pow(CurrentBeam->L0,2)
        );
    }

    // Check if the beam fractured
    if(CurrentBeam->fracture)
        if(CurrentBeam->isFractured(force, torque1)){
            // Delete the connection of the BeamConnectedSpheres. This will delete the Beam from both spheres.
            YADE_PTR_CAST<BeamConnectedSphere>(Body::byId(id2)->shape)->delConnection(Beam_id);     
            // TO DO: Delete the Beam from the scene       
            return false;
        }

    // Add the forces and torques to the bodies. Remember to transform the force and torque back to the global frame
    force = R*force;
    scene->forces.addForce(id1, force);
    scene->forces.addTorque(id1, R*torque1);
    scene->forces.addForce(id2, -force); // Newton's third law
    scene->forces.addTorque(id2, R*torque2);

	return true;
}
} // namespace yade