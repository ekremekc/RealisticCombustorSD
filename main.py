from mpi4py import MPI
from petsc4py import PETSc
from helmholtz_x.helmholtz_pkgx.active_flame_x import ActiveFlame
from helmholtz_x.helmholtz_pkgx.flame_transfer_function_x import state_space, n_tau
from helmholtz_x.helmholtz_pkgx.eigensolvers_x import fixed_point_iteration_eps, newton_solver
from helmholtz_x.helmholtz_pkgx.passive_flame_x import PassiveFlame
from helmholtz_x.helmholtz_pkgx.eigenvectors_x import normalize_eigenvector, normalize_adjoint
from helmholtz_x.geometry_pkgx.xdmf_utils import XDMFReader
from dolfinx.io import XDMFFile

import datetime
start_time = datetime.datetime.now()

import params

# Read mesh 

combustor = XDMFReader("MeshDir/combustor")
mesh, subdomains, facet_tags = combustor.getAll()

t_imap = mesh.topology.index_map(mesh.topology.dim)
num_cells = t_imap.size_local + t_imap.num_ghosts
cell_number_gathered = MPI.COMM_WORLD.allreduce(num_cells, op=MPI.SUM)
if MPI.COMM_WORLD.rank == 0:
    print("Number of cores: ", MPI.COMM_WORLD.size)
    print("Number of cells: ", cell_number_gathered)
    
# Read mesh 
# mesh, subdomains, facet_tags = read_from_msh("MeshDir/Micca.msh", cell_data=True, facet_data=True, gdim=3)

FTF = n_tau(params.N3, params.tau)
# ________________________________________________________________________________
# EVERYWHERE İS NEUMANN EXCEPT OUTLET
#6: 'Neumann',#burner back might cause the problem

boundary_conditions = {1: 'Neumann',
                       2: 'Neumann',
                       3: 'Neumann',
                       4: 'Neumann',
                       5: 'Neumann',
                       7: 'Neumann',
                       8: 'Neumann',
                       9: 'Neumann',
                       10: 'Neumann',
                       11: 'Neumann',
                       12: 'Neumann',
                       13: 'Neumann',
                       14: 'Dirichlet',
                       15: 'Dirichlet',
                       16: 'Dirichlet'}

degree = 2

target_dir = PETSc.ScalarType(2500+30)
target_adj = PETSc.ScalarType(2500-30)
c = params.c
if MPI.COMM_WORLD.rank == 0:
    print("Sound speed has imported.")
matrices = PassiveFlame(mesh, facet_tags, boundary_conditions,
                        c=c,
                        degree=degree)

matrices.assemble_A()
matrices.assemble_C()
# A = matrices.A
# C = matrices.C
if MPI.COMM_WORLD.rank == 0:
    print("Passive matrices done.")
D = ActiveFlame(mesh, subdomains, params.x_r, params.rho_amb, params.Q_tot, params.U_bulk, FTF, degree=degree)

D.assemble_submatrices('direct')

if MPI.COMM_WORLD.rank == 0:
    print("All flames are assembled. Eigenvalue solver is starting..\n")

E = fixed_point_iteration_eps(matrices, D, target_dir**2, i=0, tol=1e-1)

if MPI.COMM_WORLD.rank == 0:
    print("Eigenvalues are found:\n")

omega_1, p_1 = normalize_eigenvector(mesh, E, i=0, degree=degree)
omega_2, p_2 = normalize_eigenvector(mesh, E, i=1, degree=degree)

if MPI.COMM_WORLD.rank == 0:
    print("Direct Eigenvalues -> ", omega_1," =? ", omega_2)

# Save eigenvectors

if MPI.COMM_WORLD.rank == 0:
    print("Direct Eigenvectors are saving now ... \n")

p_1.name = "P_1_Direct"
p_2.name = "P_2_Direct"

with XDMFFile(MPI.COMM_WORLD, "Results/p_1.xdmf", "w", encoding=XDMFFile.Encoding.HDF5) as xdmf:
    xdmf.write_mesh(mesh)
    xdmf.write_function(p_1)
with XDMFFile(MPI.COMM_WORLD, "Results/p_2.xdmf", "w", encoding=XDMFFile.Encoding.HDF5) as xdmf:
    xdmf.write_mesh(mesh)
    xdmf.write_function(p_2)
# ________________________________________________________________________________

D.assemble_submatrices('adjoint')

E_adj = fixed_point_iteration_eps(matrices, D, target_adj**2, i=0, tol=1e-1, problem_type='adjoint',print_results=False)

omega_adj_1, p_adj_1 = normalize_eigenvector(mesh, E_adj, i=0, degree=degree)
omega_adj_2, p_adj_2 = normalize_eigenvector(mesh, E_adj, i=1, degree=degree)

if MPI.COMM_WORLD.rank == 0:
    print("Adjoint Eigenvalues -> ", omega_adj_1," =? ", omega_adj_2)

p_adj_norm_1 = normalize_adjoint(omega_1, p_1, p_adj_1, matrices, D)
p_adj_norm_2 = normalize_adjoint(omega_2, p_2, p_adj_2, matrices, D)

# Save eigenvectors

p_adj_1.name = "P_1_Adjoint"
p_adj_2.name = "P_2_Adjoint"

with XDMFFile(MPI.COMM_WORLD, "Results/p_adj_1.xdmf", "w", encoding=XDMFFile.Encoding.HDF5) as xdmf:
    xdmf.write_mesh(mesh)
    xdmf.write_function(p_adj_1)
with XDMFFile(MPI.COMM_WORLD, "Results/p_adj_2.xdmf", "w", encoding=XDMFFile.Encoding.HDF5) as xdmf:
    xdmf.write_mesh(mesh)
    xdmf.write_function(p_adj_2)

if MPI.COMM_WORLD.rank == 0:
    print("Total Execution Time: ", datetime.datetime.now()-start_time)


if MPI.COMM_WORLD.rank == 0:
    print("Shape Derivatives are calculating....\n")

from helmholtz_x.geometry_pkgx.shape_derivatives_x import ShapeDerivativesDegenerate

omega = (omega_1 + omega_2)/2

shape_derivatives = ShapeDerivativesDegenerate(combustor, boundary_conditions, omega, 
                                p_1, p_2, p_adj_norm_1, p_adj_norm_2, c)
if MPI.COMM_WORLD.rank == 0:
	print(shape_derivatives)

from helmholtz_x.io_utils import dict_writer

filename = "Results/shape_derivatives"
dict_writer(filename,shape_derivatives)

