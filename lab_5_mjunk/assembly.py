import scipy.sparse as sp
from scipy.sparse.linalg import spsolve
from scipy.sparse.linalg import factorized 
import plottools as pt
import numpy as np
from scipy.sparse import identity

def compute_hat_gradients(tri):
    # Compute area
    N0, N1, N2 = tri
    area=abs(0.5*np.cross(N1-N0, N2-N0))
    
    # Compute b = (1,1,1) x (x_2^1,x_2^2,x_2^3). c is similar
    ones = np.ones(3)
    b = np.cross(tri[:,1], ones)/(2*area)
    c = np.cross(ones, tri[:,0])/(2*area)

    return (area, b, c)

def assemble_mass_matrix(P, T):
    # Define constant part M_K here
    M_ref = 1/12*np.array([[2, 1, 1], 
                           [1, 2, 1], 
                           [1, 1, 2]])
    
    # Deduce number of unkowns from dimensions/shape of P
    n_p = P.shape[0]
    # Deduce number of elements from dimensions of T
    n_t = T.shape[0]
    
    # Create sparse matrix M
    M = sp.dok_matrix((n_p, n_p))
    
    for  K in range(n_t):
        # Get local to global map from T
        l2g = T[K]
        # Get triangle nodes from P
        tri = P[T[K]]
        # Unpack nodes into N1,N2,N3
        N0,N1,N2 = tri 
        # Compute area of K. 
        # Convince yourself that the following line computes the area |K|
        # area = 0.5 * abs((N1[0] - N0[0])*(N2[1] - N0[1]) -
        #                (N1[1] - N0[1])*(N2[0] - N0[0]))
        area = abs(0.5*np.cross(N1-N0, N2-N0))
        # Use area and M_ref to compute M_K
        M_K = area*M_ref
        # Instead of 2 loops we can slice out the blocks which
        # corresponds to the entries in l2g by using the funny
        # function ix_ in numpy
        M[np.ix_(l2g, l2g)] += M_K
    
    return M

def assemble_stiffness_matrix(P, T, a):
	# Create matrix as before
	# As in the mass matrix assembly:
	# Deduce number of unkowns from dimensions of P
	# number of elements from dimensions of T and sparse matrix A
	n_p = P.shape[0]
	n_t = T.shape[0]
	A = np.zeros((n_p, n_p))
	
	for  K in range(n_t):
		# Get local to global map
		loc2glo = T[K,:]
		# Get triangle coordinates and compute area
		tri = P[T[K,:]]
		N0 = tri[0]
		N1 = tri[1]
		N2 = tri[2]
		
		# Compute abar = a((N0 + N1 + N2)/3) to approximate \int_K a(x)dx
		abar = a((N0 + N1 + N2)/3)
		
		# Compute the area and the coefficient for the hat gradients
		area, b, c = compute_hat_gradients(tri)
	
		# Numpy Arrays does not behave exactly like n x 1 matrices
		# To compute the outer product b*b.T or  b*b' in Matlab notation
		# we need the np.outer function
		A_K = abar*(np.outer(b, b) + np.outer(c, c))*area
		
		# Add local element matrix to global matrix as before
		A[np.ix_(loc2glo, loc2glo)] += A_K
	return A


def assemble_load_vector(P, T, f, qr = "midpoint_2d"):
    """ Assembles the load vector """
    
    # Deduce number of unkowns from dimensions/shape of P
    n_p = P.shape[0]
    # Deduce number of elements from dimensions of T
    n_t = T.shape[0]
    
    
    # Create and intialize vector
    b = np.zeros(n_p)
    
    # Iterate over all triangles
    for  K in range(n_t):
        l2g = T[K]   # Get local to global map
        tri = P[T[K]]  # Get triangle coordinates and compute area
        N0,N1,N2 = tri 
        area = 0.5 * abs((N1[0] - N0[0])*(N2[1] - N0[1]) -
                        (N1[1] - N0[1])*(N2[0] - N0[0]))

        if qr == "midpoint_2d":   
            # 2d midpoint
            # three midpoint coordinates
            N01 = (N0 + N1)/2
            N12 = (N1 + N2)/2
            N20 = (N0 + N2)/2
            
            b_K = area/6*np.array([f(N01)+f(N20), f(N01)+f(N12), f(N12)+f(N20)])
        else:
            # 2d Trapezoid
            b_K = area/3*np.array([f(N0[0],N0[1]), f(N1[0],N1[1]), f(N2[0],N2[1])])
        # Add local contributions to the global load vector
        b[l2g] += b_K
        
    return b

def assemble_load_vector_neu(P, T, f, E, g_N, qr = "midpoint_2d"):
    """ Assembles the load vector """
    
    # Deduce number of unkowns from dimensions/shape of P
    n_p = P.shape[0]
    # Deduce number of elements from dimensions of T
    n_t = T.shape[0]
    #Deduce number of neumann edges
    n_Ne = E.shape[0]
    
    # Create and intialize vector
    b = np.zeros(n_p)
    
    # Iterate over all triangles
    for  K in range(n_t):
        l2g = T[K]   # Get local to global map
        tri = P[T[K]]  # Get triangle coordinates and compute area
        N0,N1,N2 = tri 
        area = 0.5 * abs((N1[0] - N0[0])*(N2[1] - N0[1]) -
                        (N1[1] - N0[1])*(N2[0] - N0[0]))

        if qr == "midpoint_2d":   
            # 2d midpoint
            # three midpoint coordinates
            N01 = (N0 + N1)/2
            N12 = (N1 + N2)/2
            N20 = (N0 + N2)/2
            
            b_K = area/6*np.array([f(N01)+f(N20), f(N01)+f(N12), f(N12)+f(N20)])
        else:
            # 2d Trapezoid
            b_K = area/3*np.array([f(N0[0],N0[1]), f(N1[0],N1[1]), f(N2[0],N2[1])])
        # Add local contributions to the global load vector
        b[l2g] += b_K

    for e in range(n_Ne):
        l2g = E[e]
        N0 = P[l2g[0]] #get node coordinates
        N1 = P[l2g[1]]
        edge_len = np.sqrt((N0[0]-N1[0])**2+(N0[1]-N1[1])**2)
        b_Ne = 0.5*edge_len*np.array([g_N(N0[0],N0[1]), g_N(N1[0],N1[1])])
        b[l2g] += b_Ne 
    return b



def apply_bcs_to_A(A, dirichlet_nodes):
	# Incorporate boundary conditions in matrix A
	# Set all rows corresponding to Dirichlet nodes to 0
	A[dirichlet_nodes, :] = 0 
	
	# Set diagonal to one
	A[dirichlet_nodes, dirichlet_nodes] = 1
	return A

def apply_bcs_to_b(b, g_D_values, dirichlet_nodes):
	# Incorporate boundary condition in vector b
	b[dirichlet_nodes] = g_D_values
	return b

def extract_nodes(P, inside_domain):
    return [ i for i in range(P.shape[0]) if inside_domain(P[i]) ]