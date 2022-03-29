import numpy as np
import MeshDir.dimensions as dim
from cmath import sin,cos
from dolfinx.fem import Function, FunctionSpace
from mpi4py import MPI
from helmholtz_x.geometry_pkgx.xdmf_utils import XDMFReader

def cyl2cart(rho, phi, zeta):
    # cylindrical to Cartesian
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    z = zeta
    return x, y, z


def cart2cyl(x, y, z):
    # cylindrical to Cartesian
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    zeta = z
    return rho, phi, zeta

# ------------------------------------------------------------
# flame
r_f = (dim.R_cc_inner.m+dim.R_cc_outer.m)/2 # [m]
theta = np.deg2rad(360/dim.N_s)  # [rad]
z_f = 0  # [m]

# reference
r_r = r_f
z_r = - 0.1  # [m]

# Flame transfer function
Q_tot = 800000  # [W] **per burner**
U_bulk = 18  # [m/s]

# n = Q_tot/U_bulk  # [J/m]
N3 = 1  # [/]
tau = 0.003  # [s]
rho_amb = 1.22 # [kg/m3]
# ------------------------------------------------------------
x_f = np.array([cyl2cart(r_f, i*theta, z_f) for i in range(dim.N_s)])
x_r = np.array([cyl2cart(r_r, i*theta, z_r) for i in range(dim.N_s)])

combustor = XDMFReader("MeshDir/combustor")
mesh, subdomains, facet_tags = combustor.getAll()

T_gas = 300 #K
T_amb = 300
T_end = 900
T_peak = 1300

l_cc = 0.24
r_plenum_extended_inner = 0.190
r_inner_plenum_outer = 0.120
r_cc_outer = 0.191 

alpha_theta = 360
alpha_r = 300
beta_r = 0.0198 # 0.019
gamma_theta = 1
gamma_r = 1
N_s = 20

def temperature(mesh):
    V = FunctionSpace(mesh, ("DG", 0))
    temp = Function(V)
    x_tab = V.tabulate_dof_coordinates()
    for i in range(x_tab.shape[0]):
        midpoint = x_tab[i,:]
        x = midpoint[0]
        y = midpoint[1]
        z = midpoint[2]
        r = np.sqrt(x**2+y**2)

        if z>0 and r>r_plenum_extended_inner:
            value = (T_end-T_amb)/l_cc*z+T_amb

        elif z<0 or r< r_inner_plenum_outer:
            value = 300
            
        elif z>0 and z<l_cc and r<r_cc_outer:
            A = (T_end-T_peak)*(z/l_cc)**2 + T_peak 
            B = alpha_theta*cos(N_s*np.arctan(y/x))+gamma_theta
            C = alpha_r*sin(r/beta_r)+gamma_r
            value = A+B+C
        temp.vector.setValueLocal(i, value)

    temp.x.scatter_forward()

    return temp

def sound_speed(mesh):
    temp = temperature(mesh)
    V = FunctionSpace(mesh, ("DG", 0))
    c = Function(V)
    c.x.array[:] =  20.05 * np.sqrt(temp.x.array)
    c.x.scatter_forward()
    return c


temp = temperature(mesh)
c = sound_speed(mesh)

from dolfinx.io import XDMFFile
with XDMFFile(MPI.COMM_WORLD, "data/T.xdmf", "w") as xdmf:
    xdmf.write_mesh(mesh)
    xdmf.write_function(temp)
with XDMFFile(MPI.COMM_WORLD, "data/c.xdmf", "w") as xdmf:
    xdmf.write_mesh(mesh)
    xdmf.write_function(c)
