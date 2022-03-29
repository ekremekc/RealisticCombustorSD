import os
import os_utils
import params

from datetime import datetime
start = datetime.now()


from helmholtz_solver.helmholtz_pkg.flame_transfer_function import n_tau
from helmholtz_solver.helmholtz_pkg.mshr import MeshXDMF
from helmholtz_solver.helmholtz_pkg.passive_flame import PassiveFlame
from helmholtz_solver.helmholtz_pkg.active_flame import ActiveFlame
from helmholtz_solver.helmholtz_pkg.eigensolvers import fixed_point_iteration_eps
from helmholtz_solver.helmholtz_pkg.eigenvectors import normalize_eigenvector, normalize_adjoint



(my_dir,
 mesh_dir,
 results_dir,
 eigenvectors_dir,
 pickle_dir) = os_utils.create_dirs()

msh_params = {'R_in_p': .14,
              'R_out_p': .21,
              'l_p': .07,
              'h_b': .0165,
              'l_b': .014,
              'h_pp': .00945,
              'l_pp': .006,
              'h_f': .018,
              'l_f': .006,
              'R_in_cc': .15,
              'R_out_cc': .2,
              'l_cc': .2,
              'l_ec': 0.041,
              'lc_1': 1e-1,
              'lc_2': 2e-2
              }

# EVERYWHERE İS NEUMANN EXCEPT OUTLET(COMBUSTION CHAMBER OUTLET)
boundary_conditions = {1: 'Neumann',
                       2: 'Neumann',
                       #3: 'Neumann',
                       #4: 'Neumann',
                       #5: 'Neumann',
                       #6: 'Neumann',
                       #7: 'Neumann',
                       8: 'Neumann',
                       9: 'Neumann',
                       10: 'Neumann',
                       11: 'Dirichlet',
                       12: 'Neumann',
                       13: 'Neumann',
                       14: 'Neumann',
			15: 'Neumann',
			16: 'Neumann',
			17: 'Neumann',
			18: 'Neumann'}
degree = 2

FTF = n_tau(params.N3, params.tau)
#sFTF = state_space(params.S1, params.s2, params.s3, params.s4)

foo = {'pl_rear': 1,
       'pl_outer': 2,
       'pl_inner': 3,
       'pl_front': 4,
       'b_lateral': 5,
       'b_front': 6,
       'pp_lateral': 7,
       'cc_rear': 8,
       'cc_outer': 9,
       'cc_inner': 10,
       'cc_front': 11
       }

# ________________________________________________________________________________

target = 2500

(mesh_subdir,
     eigenvectors_subdir) = os_utils.create_iter_subdirs(mesh_dir, eigenvectors_dir, i=0)

mesh_filename = 'MeshDir/combustor'

geometry = MeshXDMF(mesh_filename, write_xdmf_file=True)
geometry()
mesh = geometry.mesh
subdomains = geometry.subdomains
boundaries = geometry.boundaries

# ________________________________________________________________________________

matrices = PassiveFlame(mesh, boundaries, boundary_conditions,
                        c=params.c,
                        degree=degree)
matrices.assemble_A()
matrices.assemble_C()
# A = matrices.A
# C = matrices.C

D = ActiveFlame(mesh, subdomains, params.x_r, params.rho_amb, params.Q_tot, params.U_bulk, FTF, degree=degree)

D.assemble_submatrices('direct')

E = fixed_point_iteration_eps(matrices, D, target**2, i=0, tol=1e-4)

omega_1, p_1 = normalize_eigenvector(mesh, E, i=0, degree=degree)
omega_2, p_2 = normalize_eigenvector(mesh, E, i=2, degree=degree)

# ________________________________________________________________________________

D.assemble_submatrices('adjoint')

E_adj = fixed_point_iteration_eps(matrices, D, target**2, i=1, tol=1e-4, problem_type='adjoint')

omega_adj_1, p_adj_1 = normalize_eigenvector(mesh, E_adj, i=1, degree=degree)
omega_adj_2, p_adj_2 = normalize_eigenvector(mesh, E_adj, i=3, degree=degree)

p_adj_norm_1 = normalize_adjoint(omega_1, p_1, p_adj_1, matrices, D)
p_adj_norm_2 = normalize_adjoint(omega_2, p_2, p_adj_2, matrices, D)

# Save eigenvalues, eigenvectors and shape derivatives

os_utils.save_eigenvector(eigenvectors_subdir, 'p_1', p_1)
os_utils.save_eigenvector(eigenvectors_subdir, 'p_2', p_2)

os_utils.save_eigenvector(eigenvectors_subdir, 'p_adj_1', p_adj_norm_1)
os_utils.save_eigenvector(eigenvectors_subdir, 'p_adj_2', p_adj_norm_2)

eigs = {'omega_1': omega_1,
        'omega_2': omega_2,
        'omega_adj_1': omega_adj_1,
        'omega_adj_2': omega_adj_2}

os_utils.pickle_dump(pickle_dir, 'eigs', eigs, i=0)
os_utils.save_eigs_as_text(pickle_dir, eigs, i=0)

from helmholtz_solver.helmholtz_pkg.shape_derivatives import shape_derivatives

omega = (omega_1+omega_2)/2

results = shape_derivatives(geometry, boundary_conditions,
                                omega, (p_1, p_2), (p_adj_norm_1, p_adj_norm_2), params.c,
                                local=True)

print(results)

print("Total execution time is: ", datetime.now()-start)

