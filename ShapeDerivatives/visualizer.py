from helmholtz_x.geometry_pkgx.xdmf_utils import XDMFReader, derivatives_visualizer
from helmholtz_x.io_utils import dict_loader

# Read txt file that contains shape derivatives as a string
filename = "shape_derivatives"
shape_derivatives = dict_loader(filename)

# visualize in paraview
geometry = XDMFReader("combustor")
filename = "derivatives"
derivatives_visualizer(filename, shape_derivatives, geometry)













