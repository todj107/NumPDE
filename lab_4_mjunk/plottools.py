from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.path import Path
from matplotlib.patches import PathPatch

def plot_mesh_2d(P, T, 
                 dirichlet_nodes = None, 
                 boundary_edges  = None, 
                 plot_all_nodes=True, 
                 label_nodes=True, 
                 label_triangles=True):
    X = P[:,0]
    Y = P[:,1]

    # Get a new figure
    plt.figure()
    plt.triplot(X, Y, T.copy())
    if plot_all_nodes:
        plt.plot(X, Y, "o", markersize=8)

    if dirichlet_nodes:
        plt.plot(X[dirichlet_nodes], Y[dirichlet_nodes], "o", markersize=8)
        
    if label_nodes:
        for j, p in enumerate(P):
            plt.text(p[0], p[1], j, ha='right', color="red") # label the points
            
    if label_triangles:
        for j, s in enumerate(T):
            p = P[s].mean(axis=0)
            plt.text(p[0], p[1], '#%d' % j, ha='center', color="black") # label triangles
            
    if boundary_edges is not None :
        # Extract and plot nodes
        edge_nodes = [ ]
        for e in boundary_edges:
            edge_nodes += list(e)
        
        edge_nodes = [ ]
        for e in boundary_edges:
            edge_nodes.append((X[e[0]], Y[e[0]]))
            edge_nodes.append((X[e[1]], Y[e[1]]))
        
        codes = ([Path.MOVETO] + [Path.LINETO])*int(len(edge_nodes)/2) 
        path = Path(edge_nodes, codes)
        pathpatch = PathPatch(path, facecolor='None', edgecolor='magenta', linewidth=5.0)
        ax = plt.gca()
        ax.add_patch(pathpatch)

    plt.show()


def plot2D(X, Y, Z, triangles, title=''):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d') 
    ax.plot_trisurf(X, Y, Z, triangles=triangles.copy(), cmap=cm.viridis, linewidth=0.0)
    ax.set_title(title)
    ax.set_xlabel('$x$')
    ax.set_ylabel('$y$')
    ax.set_zlabel('$z$')
    plt.show()

def plot_comparison_2D(X, Y, f1, f2, triangles, title_f1='', title_f2=''):
    fig = plt.figure(figsize=plt.figaspect(0.33))
    ax = fig.add_subplot(1, 3, 1, projection='3d')
    ax.plot_trisurf(X, Y, f1, triangles=triangles.copy(), cmap=cm.viridis, linewidth=0.0)
    ax.set_title(title_f1)
    ax.set_xlabel('$x$')
    ax.set_ylabel('$y$')
    ax.set_zlabel('$z$')

    # Plot projected function
    ax = fig.add_subplot(1, 3, 2, projection='3d')
    ax.plot_trisurf(X, Y, f2, triangles=triangles.copy(), cmap=cm.viridis, linewidth=0.0)
    ax.set_title(title_f2)
    ax.set_xlabel('$x$')
    ax.set_ylabel('$y$')
    ax.set_zlabel('$z$')

    ax = fig.add_subplot(1, 3, 3, projection='3d')
    ax.plot_trisurf(X, Y, f1-f2, triangles=triangles.copy(), cmap=cm.viridis, linewidth=0.0)
    ax.set_title('Difference')
    ax.set_xlabel('$x$')
    ax.set_ylabel('$y$')
    ax.set_zlabel('$z$')
    plt.show()

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
			b_K = area/3*np.array([f(N0), f(N1), f(N2)])
		# Add local contributions to the global load vector
		b[l2g] += b_K
		
	return b

