#*************************************************************************
#*  Carlos Andrés del Valle Urberuaga
#*  cdelv@unal.edu.co
#*  2023
#*************************************************************************/
import yade.utils as utils
from yade.wrapper import *
from yade.minieigenHP import *
import math

def beamConnectedSphere(center, radius, dynamic=None, fixed=False, wire=False, color=None, highlight=False, material=-1, mask=1):
    """
    Create a new body with a BeamConnectedSphere shape.
    :param center: The center of the sphere.
    :param radius: The radius of the sphere.
    :param dynamic: Whether the body is dynamic or not. If None, the body will be dynamic if the material is not fixed.
    :param fixed: Whether the body is fixed or not. If True, the body will be fixed regardless of the material.
    :param wire: Whether the body is drawn as a wireframe or not.
    :param color: The color of the body. If None, a random color will be chosen.
    :param highlight: Whether the body is highlighted or not.
    :param material: The material of the body. If -1, the default material will be used.
    :param mask: The mask of the body.
    :return: The new BeamConnectedSphere body.
    """
    # Create a new body, and set its shape to a beam-connected sphere with the given radius.
    b = Body()
    b.shape = BeamConnectedSphere(radius=radius, color=color if color else utils.randomColor(), wire=wire, highlight=highlight)
    # Calculate the volume and moment of inertia of the sphere.
    V = (4. / 3) * math.pi * radius**3
    geomInert = (2. / 5.) * V * radius**2
    # Set up the rest of the body's properties.
    utils._commonBodySetup(b, V, Vector3(geomInert, geomInert, geomInert), material, pos=center, dynamic=dynamic, fixed=fixed)
    b.aspherical = False
    b.mask = mask
    return b



def beam(id1, id2, mask=1, damping=0.0, fracture=False, phi=0.0, cohesion=0.0, create_interaction=True):
    """
    Create a new beam between two bodies. The beam parameters are calculated using the spheres' parameters.
    The radius, Length, and the other geometric parameters are calculated assuming that the length is the distance between the centers of the spheres and the radius of the beam cross-section is the minimum of the spheres' radii.
    Young's modulus, Shear modulus, and density are calculated as the harmonic average of the sphere's material properties.
    The Rayleigh damping coefficients are calculated using the 2 smaller eigenmodes of the beam. The eigenmodes are calculated by solving the eigenvalue problem det(K - ω²M) = 0 where K is the stiffness matrix, and M is the mass matrix. We use the analytical solution. The matrix entries depend on the previously calculated beam parameters.
    :param id1: The id of the first body. It can be an integer or a Body object.
    :param id2: The id of the second body. It can be an integer or a Body object.
    :param mask: The mask of the beam.
    :param damping: The Rayleigh damping fraction of the beam.
    :param fracture: Whether the beam can fracture or not.
    :param phi: The angle of internal friction for the Mohr-Coulomb failure criterion.
    :param cohesion: The cohesion of the beam for the Mohr-Coulomb failure criterion.
    :param create_interaction: Whether to create an interaction between the two bodies or not.
    :return: The new beam.
    """
    if not (isinstance(id1, int) and isinstance(id2, int)):
        if isinstance(id1, Body) and isinstance(id2, Body):
            id1 = id1.id
            id2 = id2.id
        else:
            raise ValueError("beam: The ids must be integers or Body objects.")

    if id1 < 0 or id2 < 0:
        raise ValueError("beam: The ids must be positive. Probably one of the bodies hasn't been appended to O.bodies.")
    
    # Get the bodies
    sph1 = O.bodies[id1]
    sph2 = O.bodies[id2]

    # Check if the bodies are BeamConnectedSphere
    if not isinstance(sph1.shape, BeamConnectedSphere) or not isinstance(sph2.shape, BeamConnectedSphere):
        raise ValueError("beam: The bodies are not BeamConnectedSphere.")

    # Create the beam
    b = Body()
    b.shape = Beam(Damping=damping, fracture=fracture, Phi=phi, Cohesion=cohesion)
    utils._commonBodySetup(b, 1.0, Vector3(1.0, 1.0, 1.0), -1, pos=sph1.state.pos)
    b.aspherical = False
    b.mask = mask
    b.bounded = False
    
    # Create the interaction
    if create_interaction:
        utils.createInteraction(id1, id2)

    # Create the connections between the beam and the spheres
    sph1.shape.addConnection(b)
    sph2.shape.addConnection(b)

    # Configure the beam
    b.shape.configureBeam(sph1, sph2)

    return b
    

def BeamConnecNeighborSpheres(buffer = 0.0, mask=1, damping=0.0, fracture=False, phi=0.0, cohesion=0.0, create_interaction=True):
    """
    Iterate over all body pairs. If they are BeamConnectedSphere, they are not members of the same clump, there
    is no beam between them, and they are overlapping or at least closer than the buffer, then create a beam between them.
    :param buffer: The buffer distance. If the distance between the spheres is less than the sum of their radii plus the buffer and the other conditions are satisfied, then a beam will be created.
    :param mask: The mask of the beam.
    :param damping: The Rayleigh damping fraction of the beam.
    :param fracture: Whether the beam can fracture or not.
    :param phi: The angle of internal friction for the Mohr-Coulomb failure criterion.
    :param cohesion: The cohesion of the beam for the Mohr-Coulomb failure criterion.
    :param create_interaction: Whether to create an interaction between the two bodies or not.
    """
    for i in range(len(O.bodies)):
        for j in range(i+1, len(O.bodies)):
            if isinstance(O.bodies[i].shape, BeamConnectedSphere) and isinstance(O.bodies[j].shape, BeamConnectedSphere): # check that they are BeamConnectedSphere
                if O.bodies[i].clumpId != O.bodies[j].clumpId or O.bodies[j].clumpId < 0: # check that they are not from the same clump
                    if not O.bodies[i].shape.isConnectedTo(O.bodies[j]): # check that they are not already connected
                        if O.bodies[i].shape.radius + O.bodies[j].shape.radius - (O.bodies[i].state.pos - O.bodies[j].state.pos).norm() + buffer > 0: # check that they are close enough
                             O.bodies.append(beam(i, j, mask=mask, damping=damping, fracture=fracture, phi=phi, cohesion=cohesion, create_interaction=create_interaction))

