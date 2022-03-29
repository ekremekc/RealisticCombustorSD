import gmsh
import os
import sys
"""
This file automated surface tagging by calculating and pointing of mass centers and mass values.
So any addition/deletion of air admission holes might break automation-e.g. combustion chamber walls.
So please check the corresponding tag manually and implement it in one of the "elif" conditions during 
looping surfaces.. 
- OUTLETS should have different tags
"""


dir_path = os.path.dirname(os.path.realpath(__file__))
print(dir_path)

def fltk_options():

    # Type of entity label (0: description,
    #                       1: elementary entity tag,
    #                       2: physical group tag)
    gmsh.option.setNumber("Geometry.LabelType", 2)

    gmsh.option.setNumber("Geometry.PointNumbers", 0)
    gmsh.option.setNumber("Geometry.LineNumbers", 0)
    gmsh.option.setNumber("Geometry.SurfaceNumbers", 2)
    gmsh.option.setNumber("Geometry.VolumeNumbers", 2)

    # Mesh coloring(0: by element type, 1: by elementary entity,
    #                                   2: by physical group,
    #                                   3: by mesh partition)
    gmsh.option.setNumber("Mesh.ColorCarousel", 0)

    gmsh.option.setNumber("Mesh.Lines", 0)
    gmsh.option.setNumber("Mesh.SurfaceEdges", 0)
    gmsh.option.setNumber("Mesh.SurfaceFaces", 0) # CHANGE THIS FLAG TO 0 TO SEE LABELS

    gmsh.option.setNumber("Mesh.VolumeEdges", 2)
    gmsh.option.setNumber("Mesh.VolumeFaces", 2)

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)

gmsh.model.add("combustor")
gmsh.option.setString('Geometry.OCCTargetUnit', 'M')

path = os.path.dirname(os.path.abspath(__file__))

filename = "combustor"

gmsh.model.occ.importShapes(os.path.join(path, 'flame.step'))
# gmsh.model.occ.synchronize()
gmsh.model.occ.importShapes(os.path.join(path, 'fusion.step'))
gmsh.model.occ.synchronize()
gmsh.model.occ.removeAllDuplicates()
gmsh.model.occ.synchronize()

flame_vol_tags=gmsh.model.getEntities(dim=3)
# print(flame_vol_tags)
for i in range(0, 20):
    gmsh.model.addPhysicalGroup(3, [flame_vol_tags[i][1]], tag=i)

gmsh.model.addPhysicalGroup(3, [flame_vol_tags[-1][1]], tag=99)

import numpy as np
import dimensions

surfaces = gmsh.model.occ.getEntities(dim=2)

plenum_inlet, plenum_inlet_mark = [], 1
plenum_outer, plenum_outer_mark = [], 2
plenum_inner, plenum_inner_mark = [], 3
plenum_back, plenum_back_mark = [], 4

burner_lateral, burner_lateral_mark = [], 5
burner_back, burner_back_mark = [], 6

injector_lateral, injector_lateral_mark = [], 7

flame_base, flame_base_mark = [], 8

air_admission_inner, air_admission_inner_mark = [] , 9
air_admission_outer, air_admission_outer_mark = [] , 10

cc_inner, cc_inner_mark = [], 11
cc_outer, cc_outer_mark = [], 12

holes, holes_mark = [], 13

inner_air_channel_outlet, inner_air_channel_outlet_mark = [],14
cc_outlet, cc_outlet_mark = [], 15
outer_air_channel_outlet, outer_air_channel_outlet_mark = [],16



mass1 = []
mass2 = []

for surface in surfaces:
    com = gmsh.model.occ.getCenterOfMass(surface[0], surface[1])
    # print(surface, com[2])
    if np.isclose(com[2], [dimensions.z_start.m]): #PLENUM INLET TAG1
        plenum_inlet.append(surface[1])

    elif np.isclose(com[2], [(dimensions.L_combustion_chamber.m+dimensions.z_start.m)/2]): 
        mass_data = gmsh.model.occ.getMass(surface[0], surface[1])
        mass1.append([mass_data, surface])

    elif np.isclose(com[2], [-(dimensions.L_burner.m+dimensions.L_pp.m)]): #PLENUM BACK FACE # TAG4
        plenum_back.append(surface[1])

    elif np.isclose(com[2], [-(dimensions.L_burner.m+2*dimensions.L_pp.m)/2]): #BURNER LATERAL SURFACE # TAG4
        burner_lateral.append(surface[1])

    elif np.isclose(com[2], [-(dimensions.L_pp.m)]): #BURNER BACK SURFACE # TAG6
        burner_back.append(surface[1])

    elif np.isclose(com[2], [-(dimensions.L_pp.m)/2]): #INJECTOR LATERAL  SURFACE # TAG7
        injector_lateral.append(surface[1])

    elif np.isclose(com[2], [(dimensions.L_combustion_chamber.m-(dimensions.L_burner.m + dimensions.L_pp.m))/2]): # air_admission_inner  SURFACE # TAG9
        air_admission_inner.append(surface[1])
    
    elif np.isclose(com[2], [0.0696118290244798]): # air_admission_outer  SURFACE # TAG 10 # ONLY THIS ONE HAS A PROBLEM WITH AUTOMATION DUE TO HOLES
        air_admission_outer.append(surface[1])

    elif np.isclose(com[2], [(dimensions.L_combustion_chamber.m)/2]): # CC_INNER  SURFACE # TAG 11 
        cc_inner.append(surface[1])
    
    elif np.isclose(com[2], [0.12053809309687458]): # CC_OUTER  SURFACE # TAG 12 # ONLY THIS ONE HAS A PROBLEM WITH AUTOMATION DUE TO HOLES
        cc_outer.append(surface[1])

    elif np.isclose(com[2], [dimensions.z_air_hole.m]): #  HOLES # TAG 13 # 
        holes.append(surface[1])
    
    elif np.isclose(com[2], [(dimensions.L_combustion_chamber.m)]): 
        mass_data2 = gmsh.model.occ.getMass(surface[0], surface[1])
        mass2.append([mass_data2, surface])
gmsh.model.occ.synchronize()

def takeFirst(elem): # take second element for sort
    return elem[0]

mass1.sort(key=takeFirst) # sort list with key
mass2.sort(key=takeFirst) # sort list with key

face_dim = 2
gmsh.model.addPhysicalGroup(face_dim, plenum_inlet, plenum_inlet_mark) # PLENUM INLET TAG 1  
gmsh.model.addPhysicalGroup(face_dim, [mass1[-1][1][1]], plenum_outer_mark) # PLENUM OUTER TAG 2  
gmsh.model.addPhysicalGroup(face_dim, [mass1[-2][1][1]], plenum_inner_mark) # PLENUM INNER TAG 3

gmsh.model.addPhysicalGroup(face_dim, plenum_back, plenum_back_mark)              # TAG 4
gmsh.model.addPhysicalGroup(face_dim, burner_lateral, burner_lateral_mark)        # TAG 5 
gmsh.model.addPhysicalGroup(face_dim, burner_back, burner_back_mark)              # TAG 6
gmsh.model.addPhysicalGroup(face_dim, injector_lateral, injector_lateral_mark)    # TAG 7
gmsh.model.addPhysicalGroup(face_dim, air_admission_inner, air_admission_inner_mark)    # TAG 9
gmsh.model.addPhysicalGroup(face_dim, air_admission_outer, air_admission_outer_mark)    # TAG 10
gmsh.model.addPhysicalGroup(face_dim, cc_inner, cc_inner_mark)    # TAG 11
gmsh.model.addPhysicalGroup(face_dim, cc_outer, cc_outer_mark)    # TAG 12
gmsh.model.addPhysicalGroup(face_dim, holes, holes_mark)    # TAG 13


gmsh.model.addPhysicalGroup(face_dim, [mass2[0][1][1]], inner_air_channel_outlet_mark)    # TAG 14
gmsh.model.addPhysicalGroup(face_dim, [mass2[2][1][1]], cc_outlet_mark)    # TAG 15
gmsh.model.addPhysicalGroup(face_dim, [mass2[1][1][1]], outer_air_channel_outlet_mark)    # TAG 16

flame_base_indices = np.arange(0,2*dimensions.N_s+1)
for ind in flame_base_indices:
    flame_base.append(mass1[ind][1][1])
gmsh.model.addPhysicalGroup(face_dim, flame_base, flame_base_mark) # BOTTOM OF FLAME TAG 8

gmsh.option.setNumber("Mesh.MeshSizeMin", 0.03)
gmsh.option.setNumber("Mesh.MeshSizeMax", 0.03)
gmsh.model.mesh.generate(3)

gmsh.write("{}.msh".format(dir_path +"/"+filename))

if '-nopopup' not in sys.argv:
    fltk_options()
    gmsh.fltk.run()

gmsh.finalize()

from helmholtz_x.geometry_pkgx.xdmf_utils import  write_xdmf_mesh

write_xdmf_mesh(dir_path +"/"+filename,dimension=3)
