from math import *
import numpy as np
import dolfin as dolf
from scipy.io import loadmat


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

r_p = .14  # [m]
R_p = .07
l_p = .07

l_b = .014
l_pp = .006
l_f = .006
h_b = .0165
h_pp = .00945
h_f = .018
d_1 = .025
d_2 = .035

r_cc = .15
R_cc = .05
l_cc = .24

l_ec = .041  # end correction

# flame
r_f = 0.16 # [m]
theta = np.deg2rad(18)  # [rad]
z_f = 0  # [m]

# reference
r_r = r_f
z_r = - 0.08  # [m]

# ------------------------------------------------------------

r = 287.  # [J/kg/K]
gamma = 1.4  # [/]

p_amb = 101325.  # [Pa]

T_amb = 300.  # [K]

rho_amb = p_amb/(r*T_amb)  # [kg/m^3]

c_amb = sqrt(gamma*p_amb/rho_amb)  # [m/s]

T_a = 1521.  # [K] at z = 0
T_b = 1200.  # [K] at z = l_cc

# ------------------------------------------------------------

# Flame transfer function

Q_tot = 500000  # [W] **per burner**
U_bulk = 30  # [m/s]

# n = Q_tot/U_bulk  # [J/m]
N3 = 1  # [/]

tau = 0.003  # [s]



# ------------------------------------------------------------

# x_f = np.array([cyl2cart(r_f, 0*theta, z_f)])
x_f = np.array([cyl2cart(r_f, i*theta, z_f) for i in range(20)])

# x_r = np.array([cyl2cart(r_r, 0*theta, z_r)])
x_r = np.array([cyl2cart(r_r, i*theta, z_r) for i in range(20)])

# ------------------------------------------------------------
from helmholtz_solver.helmholtz_pkg.mshr import MeshXDMF

mesh_filename = 'MeshDir/combustor'

geometry = MeshXDMF(mesh_filename, write_xdmf_file=True)
geometry()
mesh = geometry.mesh

l_cc = 0.24

T_amb = 300
T_peak = 1300
T_end = 900

# radial direction 1.parantesis in expression

# azimuthal direction 2.parantesis in expression
r = 350 # alpha_theta
cons = 1 # gamma_theta
N_s = 20
# radial direction 3.parantesis in expression
#
radial_amplitude = 300 # alpha_r
period = 0.019 # beta_r
increment = 1 # gamma_r

T = dolf.Expression('''
	 x[2] >=0 &&  sqrt(pow(x[0],2)+pow(x[1],2))>0.190  ?  (T_end-T_amb)/l_cc*x[2] + T_amb :
         x[2] <= 0 || sqrt(pow(x[0],2)+pow(x[1],2))<0.130 ? 300. :
         x[2] <= l_cc && sqrt(pow(x[0],2)+pow(x[1],2))<0.191 ? ((T_end - T_peak) * pow(x[2]/l_cc, 2) + T_peak)+(r*sin(N_s*(atan(x[1]/x[0]))+cons)) +(radial_amplitude*sin(sqrt(pow(x[0],2)+pow(x[1],2))/period)+increment):
                    T_end
                    ''',
                    degree=2, l_cc=l_cc, T_peak = T_peak, T_end = T_end,T_amb=T_amb, r=r,cons = cons, radial_amplitude=radial_amplitude,period=period,increment=increment,N_s=N_s)


c = dolf.Expression('''20.05*sqrt(T)''',degree=2, T=T)
          
V = dolf.FunctionSpace(mesh,'CG',1)
u = dolf.Function(V)
u2 = dolf.Function(V)
                    
u.interpolate(T)
u2.interpolate(c)

file_to_save   = dolf.File("data/Temperature.pvd")
file_to_save   << u      

file_to_save2  = dolf.File("data/c.pvd")
file_to_save2  << u2              
