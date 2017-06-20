import numpy as np
import scipy.linalg
import scipy.sparse
import scipy.sparse.linalg
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt
import matplotlib
import scipy.interpolate
import netCDF4

import warnings
warnings.simplefilter(action = "ignore", category = FutureWarning)


class TwoDimCoupledFEM():
    
    def __init__(self, nodes, connect, mu, nu, alpha, bc_nodes=[], bc_vals=[]):
        
        self.X = nodes[:,0]
        self.Y = nodes[:,1]
        self.connect = (connect - 1)
        
        self.nu = nu
        self.mu = mu
        self.alpha = alpha
        
        self.num_elem = len(self.connect)
        
        
    def dNdxi(self, xi, eta):
        
        return [eta/4. - eta**2/4. - (eta*xi)/2. + (eta**2*xi)/2.,
                -eta/4. + eta**2/4. - (eta*xi)/2. + (eta**2*xi)/2.,
                eta/4. + eta**2/4. + (eta*xi)/2. + (eta**2*xi)/2.,
                -eta/4. - eta**2/4. + (eta*xi)/2. + (eta**2*xi)/2.,
                eta*xi - eta**2*xi,
                0.5 - eta**2/2. + xi - eta**2*xi,
                -(eta*xi) - eta**2*xi,
                -0.5 + eta**2/2. + xi - eta**2*xi,
                -2*xi + 2*eta**2*xi]
                
                
    def dNdeta(self, xi, eta):
        
        return [xi/4. - (eta*xi)/2. - xi**2/4. + (eta*xi**2)/2.,
                -xi/4. + (eta*xi)/2. - xi**2/4. + (eta*xi**2)/2.,
                xi/4. + (eta*xi)/2. + xi**2/4. + (eta*xi**2)/2.,
                -xi/4. - (eta*xi)/2. + xi**2/4. + (eta*xi**2)/2.,
                -0.5 + eta + xi**2/2. - eta*xi**2,
                -(eta*xi) - eta*xi**2,
                0.5 + eta - xi**2/2. - eta*xi**2,
                eta*xi - eta*xi**2,
                -2*eta + 2*eta*xi**2]
    
    def Nu(self, xi, eta):
                
        return [(eta*xi)/4. - (eta**2*xi)/4. - (eta*xi**2)/4. + (eta**2*xi**2)/4.,
                -(eta*xi)/4. + (eta**2*xi)/4. - (eta*xi**2)/4. + (eta**2*xi**2)/4.,
                (eta*xi)/4. + (eta**2*xi)/4. + (eta*xi**2)/4. + (eta**2*xi**2)/4.,
                -(eta*xi)/4. - (eta**2*xi)/4. + (eta*xi**2)/4. + (eta**2*xi**2)/4.,
                -eta/2. + eta**2/2. + (eta*xi**2)/2. - (eta**2*xi**2)/2.,
                xi/2. - (eta**2*xi)/2. + xi**2/2. - (eta**2*xi**2)/2.,
                eta/2. + eta**2/2. - (eta*xi**2)/2. - (eta**2*xi**2)/2.,
                -xi/2. + (eta**2*xi)/2. + xi**2/2. - (eta**2*xi**2)/2.,
                1 - eta**2 - xi**2 + eta**2*xi**2]
    
    
    def Np(self, xi, eta):
        
        return [0.25 * (1 - xi) * (1 - eta),   # - -
                0.25 * (1 - xi) * (1 + eta),   # - +
                0.25 * (1 + xi) * (1 + eta),   # + +
                0.25 * (1 + xi) * (1 - eta) ]  # - + -
                
    def dNpdxi(self, xi, eta):
        
        return [-1*(1-eta)/4, (1-eta)/4, (1+eta)/4, -1*(1+eta)/4]

    def dNpdeta(self, xi, eta):
        
        return [-1*(1-xi)/4, -1*(1-xi)/4, (1+xi)/4, (1-xi)/4]        
        
    def dNp(self, xi, eta):

        dNp = np.zeros([2,4])
        dNp[0,:] = self.dNpdxi(xi, eta)
        dNp[1,:] = self.dNpdeta(xi, eta)
        
        return dNp                        

    def compute_tensor_mobility(self, k, mu):
        
        # Mobility tensor is 2x2 (2D example)
        
        Q = np.zeros([2,2])
        Q[0,0] = k/mu
        Q[1,1] = k/mu        
        
        return Q
    
    def compute_jacobian_matrix_and_inverse(self, xi, eta):
        
        x = self.X
        y = self.Y
        con = self.connect
        
        J11 = np.dot(x[con], self.dNdxi(xi, eta))
        J12 = np.dot(y[con], self.dNdxi(xi, eta))
        J21 = np.dot(x[con], self.dNdeta(xi, eta))
        J22 = np.dot(y[con], self.dNdeta(xi, eta))
        
        self.detJ = J11 * J22 - J12 * J21
        
        self.Jinv11 =  J22 / self.detJ
        self.Jinv12 = -J12 / self.detJ
        self.Jinv21 = -J21 / self.detJ
        self.Jinv22 =  J11 / self.detJ
        
        
    def compute_B_matrix(self, xi, eta):
        
        self.compute_jacobian_matrix_and_inverse(xi, eta)
        
        dNxi = self.dNdxi(xi, eta)
        dNeta = self.dNdeta(xi, eta)
        
        Nmat = np.zeros((4, 2 * len(dNxi)))
        Nmat[0,0::2] = dNxi
        Nmat[1,0::2] = dNeta
        Nmat[2,1::2] = dNxi
        Nmat[3,1::2] = dNeta
        
        zero = np.zeros(len(self.detJ))
        
        Jmat = np.array([[self.Jinv11, self.Jinv12, zero, zero],
                         [self.Jinv21, self.Jinv22, zero, zero],
                         [zero, zero, self.Jinv11, self.Jinv12],
                         [zero, zero, self.Jinv21, self.Jinv22]])
        
        Dmat = np.array([[1.,0.,0.,0.],[0.,0.,0.,1.],[0.,1.,1.,0.]])
        
        #B = D * J * N
        return np.einsum('ij,jk...,kl',Dmat,Jmat,Nmat)
        
        
    def compute_stiffness_integrand(self, xi, eta, nu, mu):
            
        Ey = 2 * mu * (1 + nu)
            
        c11 = Ey * (1 - nu * nu) / ((1 + nu) * (1 - nu - 2 * nu * nu))
        c12 = Ey * nu / (1 - nu - 2 * nu * nu)
        c66 = Ey / (2 * (1 + nu))
        
        # Strain/Displacement Matrix             
        Cmat = np.array([[c11, c12, 0], [c12, c11, 0], [0, 0, c66]]);
            
        self.Bmat = self.compute_B_matrix(xi, eta)
            
        #K_{il} = B_{ji} C_{jk} B_{kl} \det(J)
        return np.einsum('...ji,jk,...kl,...',self.Bmat,Cmat,self.Bmat,self.detJ)
    
    
    def compute_coupling_integrand(self, xi, eta, alpha):
        
        #Q = B^T * m *alpha * Np
        m = np.array([1.0, 1.0, 0.0]) * alpha
        
        Np = self.Np(xi, eta)
        
        return np.einsum('...ij,i,k,...',self.Bmat,m,Np,self.detJ)

    def compute_storativity_integrand(self, xi, eta, alpha, phi, Ks, Kf):
    
        #S = Np^T * ( (alpha - phi)/Ks + phi/Kf ) * Np
        
        Np = self.Np(xi, eta)
        stor = (alpha-phi)/Ks + phi/Kf
        
        return np.einsum('i,...,k,...',Np,stor,Np,self.detJ)
   
   
    def compute_permeability_integrand(self, xi, eta, Q):
        
        # H = dNp^T * (k/mu) * dNp
        
        dNp = self.dNp(xi, eta)
        
        return np.einsum('...ji,jk,...kl,...',dNp,Q,dNp,self.detJ)
            
    def compute_body_force_RHS(self, xi, eta, rho, g):
    
        # Fg = Nu^T * rho * g 
        
        Nu = self.Nu(xi, eta)
        Nu = np.array(Nu).reshape([1,np.size(Nu)])
        
        return np.einsum('...ji,...,...,...',Nu,rho,g,self.detJ)
        
    def compute_fluid_pressure_force(self, xi, eta, Q, rho_f):
    
        # Fp = dNp^T * (k/mu) * rho_f * g
        g = np.zeros([2,1])
        g[1] = self.g
        dNp = self.dNp(xi, eta)
        return np.einsum('...ji,jk,...,kl...,...',dNp,Q,rho_f,g,self.detJ)        
        return dNp
        
    
    def integrate_element_matrices(self):
            
        #Use 3 x 3 Gauss integration
        wts = [5 / 9., 8 / 9., 5 / 9.]
        pts = [-np.sqrt(3 / 5.), 0.0, np.sqrt(3 / 5.)]
        
        # Make up permeability/viscosity info
        perm    = 10. *9.8692327E-16            # m**2
        visc    = 0.002                         # Pa*sec
        phi     = 0.1
        Vp      = 2360                          # m/s
        Vs      = 1302                          # m/s
        rho_m   = 2260                          # kg/m**3
        rho_f   = 1000                          # kg/m**3
        rho     = rho_m*(1.0-phi) + rho_f*phi   # Bulk Density, kg/m**3
        Ks      = rho*(Vp**2 - 0.75 * Vs**2)    # Bulk (solid) Modulus, Pa
        Kf      = 1/1e9                         # Fluid Modulus, Pa (1/cf)
        
        
        
        # Ditto Delta T
        
        self.delta_t = 3600         # sec       
        
        
        mob_tensor = self.compute_tensor_mobility(perm,visc)
        # Temporary global setting
        self.mob_tensor = self.compute_tensor_mobility(perm,visc)        
        
            
        K = np.zeros((self.num_elem, 18, 18))
        Q = np.zeros((self.num_elem, 18,  4))
        H = np.zeros((self.num_elem,  4,  4))
        S = np.zeros((self.num_elem,  4,  4))
            
        for i in range(3):
            for j in range(3):
                K += wts[i] * wts[j] * self.compute_stiffness_integrand(pts[i], pts[j], self.nu, self.mu)
                Q += wts[i] * wts[j] * self.compute_coupling_integrand(pts[i], pts[j], self.alpha)
                H += wts[i] * wts[j] * self.compute_permeability_integrand(pts[i], pts[j], mob_tensor)
                S += wts[i] * wts[j] * self.compute_storativity_integrand(pts[i], pts[j], self.alpha, phi, Ks, Kf)
                                    
        return (K,Q,H,S)
    
    def integrate_RHS(self):

        #Use 3 x 3 Gauss integration
        wts = [5 / 9., 8 / 9., 5 / 9.]
        pts = [-np.sqrt(3 / 5.), 0.0, np.sqrt(3 / 5.)]
        
        Fgu = np.zeros((self.num_elem,9,1))
        Fgp = np.zeros((self.num_elem,4,1))    
        
        for i in range(3):
            for j in range(3):
                Fgu += wts[i] * self.compute_body_force_RHS(pts[i],pts[j],self.rho,self.g)
                Fgp += wts[i] * self.compute_fluid_pressure_force(pts[i], pts[j], self.mob_tensor, self.rho_f)
        return Fgu, Fgp   

    
    def set_mesh_properties(self):

        # Make up permeability/viscosity info
        self.perm    = 10. *9.8692327E-16            # m**2
        self.visc    = 0.002                         # Pa*sec
        self.phi     = 0.1
        self.Vp      = 2360                          # m/s
        self.Vs      = 1302                          # m/s
        self.rho_m   = 2260                          # kg/m**3
        self.rho_f   = 1000                          # kg/m**3
        self.rho     = self.rho_m*(1.0-self.phi) + self.rho_f*self.phi   # Bulk Density, kg/m**3
        self.Ks      = self.rho*(self.Vp**2 - 0.75 * self.Vs**2)    # Bulk (solid) Modulus, Pa
        self.Kf      = 1/1e9                         # Fluid Modulus, Pa (1/cf)
        self.g       = 9.80                          # m/s**2                
    
    def compute_stress_at_gauss_point(self, xi, eta, disp):
        
        mu = self.mu
        nu = self.nu
        
        Ey = 2 * mu * (1 + nu)
            
        c11 = Ey * (1 - nu * nu) / ((1 + nu) * (1 - nu - 2 * nu * nu))
        c12 = Ey * nu / (1 - nu - 2 * nu * nu)
        c66 = Ey / (2 * (1 + nu))
            
        Cmat = np.array([[c11, c12, 0], [c12, c11, 0], [0, 0, c66]]);
            
        Bmat = self.compute_B_matrix(xi, eta)
        
        elem_disp = disp[self.connect].reshape(-1,18)
        
        #stress_{i} = C_{ik} B_{kl} disp{l}
        return np.einsum('ik,...kl,...l',Cmat,Bmat,elem_disp).reshape(-1,3)
    
    
    def compute_stress(self, disp):
            
        #Gauss points
        pts = [-np.sqrt(3 / 5.), 0.0, np.sqrt(3 / 5.)]
            
        return np.array([[ self.compute_stress_at_gauss_point(i, j, disp) for i in pts ] for j in pts])
    
    
    def compute_gauss_point_locations(self, coords):
            
        #Gauss points
        pts = [-np.sqrt(3 / 5.), 0.0, np.sqrt(3 / 5.)]
        
        X = coords[:,0][self.connect]
        Y = coords[:,1][self.connect]
            
        xloc = np.array([[ np.dot(X, self.Nu(i, j)) for i in pts ] for j in pts]).flatten()
        yloc = np.array([[ np.dot(Y, self.Nu(i, j)) for i in pts ] for j in pts]).flatten()
        
        return (xloc, yloc)
    
    
    def assemble(self):
        
        # Institute debug system values
        self.set_mesh_properties()
        
        #Construct a DOF map.  We'll start by assuming that every node has 3 DOF
        #this is obviously not true, but we'll correct it.
        fake_dof_map = np.zeros(3*len(self.X), dtype=np.int64).reshape(-1, 3)
        
        #The nodes that actually have pressue DOF's are the "corner" nodes of the element
        #or the first 4 nodes in each row of the connectivity array, let's select the others
        no_pressure_dof_nodes = (self.connect[:,4:]).flatten()
        
        #Now for these nodes that do not have pressure DOF's, we'll "flag" those DOFs with a -1
        fake_dof_map[:,2][no_pressure_dof_nodes] = -1
        fake_dof_map = fake_dof_map.flatten()
        
        #Now the rest of the DOF indices are not monotonically increasing (we removed some of them
        #and replaced them with -1's), so now we need to replace the non -1 entries with a monotonically
        #increasing DOF map that corresponds to the total number of DOF's.  First let's figure out
        #how many total DOF's ther are, this corresponds to wherever there are 0's in the fake_dof_map
        total_dof = (fake_dof_map == 0).sum()
        
        #Create a monotonically increasing range from 0 to the total_dof
        dof_range = np.arange(total_dof, dtype=np.int64)
        
        #Replace the 0's in the fake_dof_map with the monotonically increasing range
        fake_dof_map[np.where(fake_dof_map != -1)] = dof_range
        
        #Create the real dof_map, there will still be -1 "flags" in the third column
        self.dof_map = fake_dof_map.reshape(-1,3)
        
        #Allocate global stiffness matrix and r.h.s vector
        self.Mat_K = np.zeros((total_dof, total_dof), dtype=np.double)
        self.Mat_K_n = np.zeros((total_dof, total_dof), dtype=np.double)        
        self.Vec_F = np.zeros(total_dof, dtype=np.double)
        
        #The DOF indices cooresponding to displacement
        id_disp = (self.dof_map[:,:2][self.connect]).reshape(-1,18)
        # DOF indicies for x dir displacement
        idx_disp = (self.dof_map[:,0][self.connect]).reshape(-1,9)                
        # DOF indicies for x dir displacement
        idy_disp = (self.dof_map[:,1][self.connect]).reshape(-1,9)                
        
        #The DOF indices cooresponding to pressure, they should not have any -1's
        #because we only choose those that in the first 4 columns of each row in
        #the connectivity
        id_pres = self.dof_map[:,-1][self.connect[:,:4]]
        
        #Integrate element stiffness matrices
        K, Q, H, S = self.integrate_element_matrices()

        #Integrate RHS
        self.Fgu_old = np.zeros((self.num_elem,9,1))
        
        Fgu, Fgp = self.integrate_RHS()        
                
        #Assemble into global stiffness matrix
        for i in range(self.num_elem):
            
            # Add elastic stiffness matrix            
            idx_grid_disp = np.ix_(id_disp[i], id_disp[i])
            self.Mat_K[idx_grid_disp]  += -1*K[i]
            
            # Add coupling matrix
            idx_grid_pres = np.ix_(id_disp[i], id_pres[i])
            self.Mat_K[idx_grid_pres] += Q[i]
            
            # Add transposed coupling matrix
            idx_grid_pres = np.ix_(id_pres[i], id_disp[i])
            self.Mat_K[idx_grid_pres] += Q[i].T

            # Add Storativity Matrix            
            idx_grid_pres = np.ix_(id_pres[i], id_pres[i])
            self.Mat_K[idx_grid_pres] += S[i]

            # Add Fluid Permeability Matrix
            idx_grid_pres = np.ix_(id_pres[i], id_pres[i])
            self.Mat_K[idx_grid_pres] += H[i] * self.delta_t / 2


        #Assemble into global stiffness matrix for RHS
            
            # Add elastic stiffness matrix            
            idx_grid_disp = np.ix_(id_disp[i], id_disp[i])
            self.Mat_K_n[idx_grid_disp]  += -1*K[i]
            
            # Add coupling matrix
            idx_grid_pres = np.ix_(id_disp[i], id_pres[i])
            self.Mat_K_n[idx_grid_pres] += Q[i]
            
            # Add transposed coupling matrix
            idx_grid_pres = np.ix_(id_pres[i], id_disp[i])
            self.Mat_K_n[idx_grid_pres] += Q[i].T

            # Add Storativity Matrix            
            idx_grid_pres = np.ix_(id_pres[i], id_pres[i])
            self.Mat_K_n[idx_grid_pres] += S[i]

            # Add Fluid Permeability Matrix
            idx_grid_pres = np.ix_(id_pres[i], id_pres[i])
            self.Mat_K_n[idx_grid_pres] -= H[i] * self.delta_t / 2

        #Institute Body Forces (grav)
        
            # Add displacement body force
            idy_vec_disp = np.ix_(idy_disp[i])
            self.Vec_F[idy_vec_disp] -= (Fgu[i].ravel() - self.Fgu_old[i].ravel()) / self.delta_t
            
            idx_vec_pres = np.ix_(id_pres[i])
            self.Vec_F[idx_vec_pres] += Fgp[i].ravel()
            
    def apply_essential_bc(self, nodes, values, dof="x"):
        
        node_idx = nodes - 1
        
        if dof == "x":
            dof_idx = 0
        elif dof == "y":
            dof_idx = 1
        elif dof == "p":
            dof_idx = 2
        
        row_replace = np.zeros(len(self.Mat_K))
        
        for value_idx, node in enumerate(node_idx):        
            
            row_idx = self.dof_map[node][dof_idx]
            
            self.Mat_K[row_idx] = row_replace
            self.Mat_K[row_idx,row_idx] = 1
            self.Vec_F[row_idx] = values[value_idx]
            
            
    def solve(self):
        
        self.Mat_K = scipy.sparse.csr_matrix(self.Mat_K)
        self.Mat_K_n = scipy.sparse.csr_matrix(self.Mat_K)
        self.Vec_F = self.Mat_K_n*self.Vec_F + self.delta_t*self.Vec_F
        
        return scipy.sparse.linalg.spsolve(self.Mat_K,self.Vec_F)
################################################################################
        
coords = np.loadtxt("coords.csv", delimiter=',', dtype=np.double)
connect = np.loadtxt("connect.csv", delimiter=',', dtype=np.int64)
ns1 = np.loadtxt("nodeset1.csv", delimiter=',', dtype=np.int64)
ns2 = np.loadtxt("nodeset2.csv", delimiter=',', dtype=np.int64)
ns3 = np.loadtxt("nodeset3.csv", delimiter=',', dtype=np.int64)
ns4 = np.loadtxt("nodeset4.csv", delimiter=',', dtype=np.int64)

# Load from .exo file

nc = netCDF4.Dataset('Inclusion2D_Coarse.exo')
connect1 =  nc.variables['connect1'][:]
connect2 =  nc.variables['connect2'][:]
connect3 =  nc.variables['connect3'][:]

connect = np.vstack( (connect1,connect2,connect3))
coordx = nc.variables['coordx'][:]
coordy = nc.variables['coordy'][:]
coords = np.vstack((coordx,coordy)).T

left    = nc.variables['node_ns2'][:]
floor   = nc.variables['node_ns3'][:]
right   = nc.variables['node_ns4'][:]
top     = nc.variables['node_ns5'][:]


problem = TwoDimCoupledFEM(coords, connect, nu=0.3, mu=9722223330304.0, alpha=0.8)

problem.assemble()

problem.apply_essential_bc(floor,np.zeros(len(floor)),dof="y")

#problem.apply_essential_bc(floor,np.ones(len(floor))*100,dof="p")
problem.apply_essential_bc(left,np.zeros(len(left)),dof="x")
problem.apply_essential_bc(right,np.zeros(len(right)),dof="x")

problem.apply_essential_bc(top,np.zeros(len(top)),dof="p")
#problem.apply_essential_bc(ns4,np.zeros(len(ns4)),dof="x")

#ns1 includes displacement nodes on the interior, we can't apply pressure
#to these, so first we modify ns1 to include only pressure nodes
ns1_pres = np.intersect1d(ns1,connect[:,:4].flatten())
#problem.apply_essential_bc(ns1_pres,5*np.ones(len(ns1_pres)),dof="p")
#Same for ns2
ns2_pres = np.intersect1d(ns2,connect[:,:4].flatten())
#problem.apply_essential_bc(ns2_pres,500*6894*np.ones(len(ns2_pres)),dof="p")

x = problem.solve()

displacement = x[problem.dof_map[:,:2].flatten()].reshape(-1,2)
deformed_pos = coords + displacement

pres_dof = problem.dof_map[:,-1][np.where(problem.dof_map[:,-1] != -1)]
pressure = x[pres_dof]


################################################################################        
# Plots

patches = []
for coord in coords[connect[:,0:4]-1]:
    quad = Polygon(coord, facecolor='none', fill=False)
    patches.append(quad)

pres_idx, = np.where(problem.dof_map[:,-1] != -1)
X = coords[pres_idx,0]
Y = coords[pres_idx,1]
grid_x, grid_y = np.mgrid[0:50:1000j, 0:-50:1000j]
p = scipy.interpolate.griddata((X, Y), pressure, (grid_x, grid_y), method='cubic')

X = coords[:,0]
Y = coords[:,1]
u_x = scipy.interpolate.griddata((X, Y), displacement[:,0], (grid_x, grid_y), method='cubic')
u_y = scipy.interpolate.griddata((X, Y), displacement[:,1], (grid_x, grid_y), method='cubic')
displacement_mag = np.sqrt(displacement[:,0] * displacement[:,0] + displacement[:,1] * displacement[:,1])
disp_mag = scipy.interpolate.griddata((X, Y), displacement_mag, (grid_x, grid_y), method='cubic')
interior = np.sqrt((grid_x**2 / 0.8**2) + (grid_y**2 / 1**2)) < 1
p[interior] = np.nan

fig, ax = plt.subplots()
cs = ax.contourf(grid_x, grid_y, p, cmap="coolwarm") #,levels=np.linspace(-2, 5, 50))
fig.colorbar(cs, ax=ax);
#colors = 100 * np.random.rand(len(patches))


p = PatchCollection(patches, match_original=True)
p.set_linewidth(0.1)
    
#p.set_array(np.array(colors))
ax.add_collection(p)
ax.set_xlim([0, 50])
ax.set_ylim([-50,0])
ax.set_aspect('equal')        
#        
# #Sigma xx

#stress = problem.compute_stress(displacement).reshape(-1,3)
#x_gauss_pt, y_gauss_pt = problem.compute_gauss_point_locations(deformed_pos)


#stress_x = scipy.interpolate.griddata((x_gauss_pt, y_gauss_pt), stress[:,0], (grid_x, grid_y), method='cubic')
#stress_y = scipy.interpolate.griddata((x_gauss_pt, y_gauss_pt), stress[:,1], (grid_x, grid_y), method='cubic')
#stress_xy = scipy.interpolate.griddata((x_gauss_pt, y_gauss_pt), stress[:,2], (grid_x, grid_y), method='cubic')
##stress_x[interior] = np.nan
##stress_y[interior] = np.nan
##stress_xy[interior] = np.nan


#plt.figure()
#plt.gca().set_aspect('equal')
#plt.contourf(grid_x, grid_y, stress_x, cmap="coolwarm",levels=np.linspace(-1, 2, 50))
#plt.colorbar();
#plt.title("$\sigma_{xx}$");


# Sigma yy

#plt.figure()
#plt.gca().set_aspect('equal')
#plt.contourf(grid_x, grid_y, stress_y, cmap="coolwarm",levels=np.linspace(-1,  2, 50))
#plt.colorbar();
#plt.title("$\sigma_{yy}$");

### Sigma xy

#plt.figure()
#plt.gca().set_aspect('equal')
#plt.contourf(grid_x, grid_y, stress_xy, cmap="coolwarm",levels=np.linspace(-1, 1, 50))
#plt.colorbar();
#plt.title("$\sigma_{xy}$");
#        
#        
## X Disp

#X = coords[:,0]
#Y = coords[:,1]
#disp_x = scipy.interpolate.griddata((X, Y), displacement[:,0], (grid_x, grid_y), method='cubic')
##disp_x[interior] = np.nan
#plt.figure()
#plt.gca().set_aspect('equal')
#plt.contourf(grid_x, grid_y, disp_x, cmap="coolwarm")
#plt.colorbar();
#plt.title("X displacement");

## Y Disp

#disp_y = scipy.interpolate.griddata((X, Y), displacement[:,1], (grid_x, grid_y), method='cubic')
##disp_y[interior] = np.nan
#plt.figure()
#plt.gca().set_aspect('equal')
#plt.contourf(grid_x, grid_y, disp_y, cmap="coolwarm")
#plt.colorbar();
#plt.title("Y displacement");

## Disp Magnitude

#displacement_mag = np.sqrt(displacement[:,0] * displacement[:,0] + displacement[:,1] * displacement[:,1])
#disp_mag = scipy.interpolate.griddata((X, Y), displacement_mag, (grid_x, grid_y), method='cubic')
##disp_mag[interior] = np.nan
#plt.figure()
#plt.gca().set_aspect('equal')
#plt.contourf(grid_x, grid_y, disp_mag, cmap="coolwarm")#,levels=np.linspace(0.0, 0.06, 50))
#plt.colorbar();
#plt.title("Displacement Magnitude");        
