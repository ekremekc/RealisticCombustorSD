import FreeCAD,FreeCADGui,os
import PartDesignGui

doc = FreeCAD.newDocument('Combustor')
path = "/home/ekrem/Dev/Acoustics-Dev/Shapes/RealCombustor/Hole/"

class Length:
	def __init__(self, length):
		self.m = length
		self. mm = length*1000
		self.str = str(self.mm)+ " mm"


L_combustion_chamber = Length(0.24)
L_pp = Length(0.02)
L_burner = Length(0.08)
L_plenum_inlet = Length(0.14)
L_total = Length(L_combustion_chamber.m + L_pp.m + L_burner.m + L_plenum_inlet.m)

L_flame = Length(0.03)

R_plenum_inner = Length(0.10)
R_plenum_outer = Length(0.22)

R_air_admission_inner = Length(0.11)
R_air_admission_outer =Length(0.20)

R_cc_inner = Length(0.12)
R_cc_outer = Length(0.19)
R_burner = Length(0.018)
R_pp = Length(0.010)
R_flame =  Length(0.020)

N_s = 20 # Number of sectors

#Cooling Hole parameters
tolerance = Length(0.005) # it is used for fully merging the combustion chamber and outer plenum
R_air_hole = Length(0.01) 
L_air_hole = Length(R_air_admission_outer.m -  R_cc_outer.m+tolerance.m) 

x_air_hole = Length(R_cc_outer.m - tolerance.m)
z_air_hole = Length(L_combustion_chamber.m*2/5)

R_mid_plane = Length((R_cc_inner.m+R_cc_outer.m)/2)

z_start = Length(-(L_pp.m + L_burner.m + L_plenum_inlet.m))
z_pp = Length(-(L_pp.m))


App.ActiveDocument.addObject("Part::Cylinder","Cylinder")
App.ActiveDocument.ActiveObject.Label = "Cylinder"
FreeCAD.getDocument('Combustor').getObject('Cylinder').Placement = App.Placement(App.Vector(0,0,z_start.mm),App.Rotation(App.Vector(0,0,1),0))
App.ActiveDocument.recompute()
### End command Part_Cylinder
# Gui.Selection.addSelection('Unnamed','Cylinder')

FreeCAD.getDocument('Combustor').getObject('Cylinder').Radius = R_plenum_outer.str
FreeCAD.getDocument('Combustor').getObject('Cylinder').Height = L_total.str

App.ActiveDocument.addObject("Part::Cylinder","Cylinder")
App.ActiveDocument.ActiveObject.Label = "Cylinder001"
FreeCAD.getDocument('Combustor').getObject('Cylinder001').Placement = App.Placement(App.Vector(0,0,z_start.mm),App.Rotation(App.Vector(0,0,1),0))
App.ActiveDocument.recompute()

FreeCAD.getDocument('Combustor').getObject('Cylinder001').Radius = R_plenum_inner.str
FreeCAD.getDocument('Combustor').getObject('Cylinder001').Height = L_total.str

### Begin command Part_Cut
App.activeDocument().addObject("Part::Cut","Cut")
App.activeDocument().Cut.Base = App.activeDocument().Cylinder
App.activeDocument().Cut.Tool = App.activeDocument().Cylinder001
Gui.activeDocument().Cylinder.Visibility=False
Gui.activeDocument().Cylinder001.Visibility=False
App.getDocument('Combustor').getObject('Cut').ViewObject.ShapeColor=getattr(App.getDocument('Combustor').getObject('Cylinder').getLinkedObject(True).ViewObject,'ShapeColor',App.getDocument('Combustor').getObject('Cut').ViewObject.ShapeColor)
App.getDocument('Combustor').getObject('Cut').ViewObject.DisplayMode=getattr(App.getDocument('Combustor').getObject('Cylinder').getLinkedObject(True).ViewObject,'DisplayMode',App.getDocument('Combustor').getObject('Cut').ViewObject.DisplayMode)
App.ActiveDocument.recompute()
Gui.SendMsgToActiveView("ViewFit")

##################################################################################################################
# Now We are substracting the inner wall region

L_emptied = Length(L_total.m - L_plenum_inlet.m)
z_burner = Length(-(L_pp.m + L_burner.m ))

App.ActiveDocument.addObject("Part::Cylinder","Cylinder")
App.ActiveDocument.ActiveObject.Label = "Cylinder002"
FreeCAD.getDocument('Combustor').getObject('Cylinder002').Placement = App.Placement(App.Vector(0,0,z_burner.mm),App.Rotation(App.Vector(0,0,1),0))
App.ActiveDocument.recompute()
### End command Part_Cylinder
# Gui.Selection.addSelection('Unnamed','Cylinder')

FreeCAD.getDocument('Combustor').getObject('Cylinder002').Radius = R_air_admission_outer.str
FreeCAD.getDocument('Combustor').getObject('Cylinder002').Height = L_emptied.str

App.ActiveDocument.addObject("Part::Cylinder","Cylinder")
App.ActiveDocument.ActiveObject.Label = "Cylinder003"
FreeCAD.getDocument('Combustor').getObject('Cylinder003').Placement = App.Placement(App.Vector(0,0,z_burner.mm),App.Rotation(App.Vector(0,0,1),0))
App.ActiveDocument.recompute()

FreeCAD.getDocument('Combustor').getObject('Cylinder003').Radius = R_air_admission_inner.str
FreeCAD.getDocument('Combustor').getObject('Cylinder003').Height = L_emptied.str

App.activeDocument().addObject("Part::Cut","Cut001")
App.activeDocument().Cut001.Base = App.activeDocument().Cylinder002
App.activeDocument().Cut001.Tool = App.activeDocument().Cylinder003
Gui.activeDocument().Cylinder002.Visibility=False
Gui.activeDocument().Cylinder003.Visibility=False
App.getDocument('Combustor').getObject('Cut001').ViewObject.ShapeColor=getattr(App.getDocument('Combustor').getObject('Cylinder002').getLinkedObject(True).ViewObject,'ShapeColor',App.getDocument('Combustor').getObject('Cut001').ViewObject.ShapeColor)
App.getDocument('Combustor').getObject('Cut001').ViewObject.DisplayMode=getattr(App.getDocument('Combustor').getObject('Cylinder002').getLinkedObject(True).ViewObject,'DisplayMode',App.getDocument('Combustor').getObject('Cut001').ViewObject.DisplayMode)

App.activeDocument().addObject("Part::Cut","Cut002")
App.activeDocument().Cut002.Base = App.activeDocument().Cut
App.activeDocument().Cut002.Tool = App.activeDocument().Cut001
Gui.activeDocument().Cut.Visibility=False
Gui.activeDocument().Cut001.Visibility=False
App.getDocument('Combustor').getObject('Cut002').ViewObject.ShapeColor=getattr(App.getDocument('Combustor').getObject('Cut').getLinkedObject(True).ViewObject,'ShapeColor',App.getDocument('Combustor').getObject('Cut002').ViewObject.ShapeColor)
App.getDocument('Combustor').getObject('Cut002').ViewObject.DisplayMode=getattr(App.getDocument('Combustor').getObject('Cut').getLinkedObject(True).ViewObject,'DisplayMode',App.getDocument('Combustor').getObject('Cut002').ViewObject.DisplayMode)
App.ActiveDocument.recompute()

##########################################################################################################
# Now we will generate burners


App.ActiveDocument.addObject("Part::Cylinder","Cylinder")
App.ActiveDocument.ActiveObject.Label = "Cylinder004"
App.ActiveDocument.recompute()
FreeCAD.getDocument('Combustor').getObject('Cylinder004').Radius = R_burner.str
FreeCAD.getDocument('Combustor').getObject('Cylinder004').Height = L_burner.str
Gui.SendMsgToActiveView("ViewFit")

FreeCAD.getDocument('Combustor').getObject('Cylinder004').Placement = App.Placement(App.Vector(R_mid_plane.mm,0,z_burner.mm),App.Rotation(App.Vector(0,0,1),0))

# copy burner

import Draft
_obj_ = Draft.make_polar_array(App.ActiveDocument.Cylinder004, number=N_s, angle=360.0, center=FreeCAD.Vector(0.0, 0.0, 0.0), use_link=True)
# Gui.Selection.addSelection('Combustor','Array')
_obj_.Fuse = False
Draft.autogroup(_obj_)
App.ActiveDocument.recompute()

##########################################################################################################
# Now we will generate injectors

z_pp = Length(-L_pp.m)
App.ActiveDocument.addObject("Part::Cylinder","Cylinder")
App.ActiveDocument.ActiveObject.Label = "Cylinder005"
App.ActiveDocument.recompute()
FreeCAD.getDocument('Combustor').getObject('Cylinder005').Radius = R_pp.str
FreeCAD.getDocument('Combustor').getObject('Cylinder005').Height = L_pp.str
Gui.SendMsgToActiveView("ViewFit")

FreeCAD.getDocument('Combustor').getObject('Cylinder005').Placement = App.Placement(App.Vector(R_mid_plane.mm,0,z_pp.mm),App.Rotation(App.Vector(0,0,1),0))

# copy burner

import Draft
_obj_ = Draft.make_polar_array(App.ActiveDocument.Cylinder005, number=N_s, angle=360.0, center=FreeCAD.Vector(0.0, 0.0, 0.0), use_link=True)
# Gui.Selection.addSelection('Combustor','Array')
_obj_.Fuse = False
Draft.autogroup(_obj_)
App.ActiveDocument.recompute()

##########################################################################################################
# Now we will generate EMPTY FLAMES

App.ActiveDocument.addObject("Part::Cylinder","Cylinder")
App.ActiveDocument.ActiveObject.Label = "Cylinder006"
App.ActiveDocument.recompute()
FreeCAD.getDocument('Combustor').getObject('Cylinder006').Radius = R_flame.str
FreeCAD.getDocument('Combustor').getObject('Cylinder006').Height = L_flame.str
Gui.SendMsgToActiveView("ViewFit")

FreeCAD.getDocument('Combustor').getObject('Cylinder006').Placement = App.Placement(App.Vector(R_mid_plane.mm,0,0),App.Rotation(App.Vector(0,0,1),0))

_obj_ = Draft.make_polar_array(App.ActiveDocument.Cylinder006, number=N_s, angle=360.0, center=FreeCAD.Vector(0.0, 0.0, 0.0), use_link=True)
_obj_.Fuse = False
Draft.autogroup(_obj_)
App.ActiveDocument.recompute()

###########################################################################################################
# Now we add combustion chamber

App.ActiveDocument.addObject("Part::Cylinder","Cylinder")
App.ActiveDocument.ActiveObject.Label = "Cylinder007"
App.ActiveDocument.recompute()

FreeCAD.getDocument('Combustor').getObject('Cylinder007').Radius = R_cc_outer.str
FreeCAD.getDocument('Combustor').getObject('Cylinder007').Height = L_combustion_chamber.str

App.ActiveDocument.addObject("Part::Cylinder","Cylinder")
App.ActiveDocument.ActiveObject.Label = "Cylinder008"
App.ActiveDocument.recompute()

FreeCAD.getDocument('Combustor').getObject('Cylinder008').Radius = R_cc_inner.str
FreeCAD.getDocument('Combustor').getObject('Cylinder008').Height = L_combustion_chamber.str

App.activeDocument().addObject("Part::Cut","Cut003")
App.activeDocument().Cut003.Base = App.activeDocument().Cylinder007
App.activeDocument().Cut003.Tool = App.activeDocument().Cylinder008
Gui.activeDocument().Cylinder007.Visibility=False
Gui.activeDocument().Cylinder008.Visibility=False
App.getDocument('Combustor').getObject('Cut003').ViewObject.ShapeColor=getattr(App.getDocument('Combustor').getObject('Cylinder007').getLinkedObject(True).ViewObject,'ShapeColor',App.getDocument('Combustor').getObject('Cut003').ViewObject.ShapeColor)
App.getDocument('Combustor').getObject('Cut003').ViewObject.DisplayMode=getattr(App.getDocument('Combustor').getObject('Cylinder007').getLinkedObject(True).ViewObject,'DisplayMode',App.getDocument('Combustor').getObject('Cut003').ViewObject.DisplayMode)
App.ActiveDocument.recompute()

############################################################################
# NOW we are exctracting flames from Combustion Chamber

App.activeDocument().addObject("Part::Cut","Cut004")
App.activeDocument().Cut004.Base = App.activeDocument().Cut003
App.activeDocument().Cut004.Tool = App.activeDocument().Array002
Gui.activeDocument().Cut003.Visibility=False
Gui.activeDocument().Array002.Visibility=False
App.getDocument('Combustor').getObject('Cut004').ViewObject.ShapeColor=getattr(App.getDocument('Combustor').getObject('Cut003').getLinkedObject(True).ViewObject,'ShapeColor',App.getDocument('Combustor').getObject('Cut004').ViewObject.ShapeColor)
App.getDocument('Combustor').getObject('Cut004').ViewObject.DisplayMode=getattr(App.getDocument('Combustor').getObject('Cut003').getLinkedObject(True).ViewObject,'DisplayMode',App.getDocument('Combustor').getObject('Cut004').ViewObject.DisplayMode)
App.ActiveDocument.recompute()

################################################################################
# Now we add cooling holes 


App.ActiveDocument.addObject("Part::Cylinder","Cylinder")
App.ActiveDocument.ActiveObject.Label = "Cylinder009"
App.ActiveDocument.recompute()
FreeCAD.getDocument('Combustor').getObject('Cylinder009').Radius = R_air_hole.str
FreeCAD.getDocument('Combustor').getObject('Cylinder009').Height = L_air_hole.str
Gui.SendMsgToActiveView("ViewFit")

FreeCAD.getDocument('Combustor').getObject('Cylinder009').Placement = App.Placement(App.Vector(x_air_hole.mm,0,z_air_hole.mm),App.Rotation(App.Vector(0,1,0),90))
_obj_ = Draft.make_polar_array(App.ActiveDocument.Cylinder009, number=N_s, angle=360.0, center=FreeCAD.Vector(0.0, 0.0, 0.0), use_link=True)
_obj_.Fuse = True
Draft.autogroup(_obj_)
App.ActiveDocument.recompute()

################################################
#Now fuse everything except flame region
App.activeDocument().addObject("Part::MultiFuse","Fusion")
App.activeDocument().Fusion.Shapes = [App.activeDocument().Cut002,App.activeDocument().Array,App.activeDocument().Array001,App.activeDocument().Cut004,App.activeDocument().Array003,]
Gui.activeDocument().Cut002.Visibility=False
Gui.activeDocument().Array.Visibility=False
Gui.activeDocument().Array001.Visibility=False
Gui.activeDocument().Cut004.Visibility=False
Gui.activeDocument().Array003.Visibility=False
App.getDocument('Combustor').getObject('Fusion').ViewObject.ShapeColor=getattr(App.getDocument('Combustor').getObject('Cut002').getLinkedObject(True).ViewObject,'ShapeColor',App.getDocument('Combustor').getObject('Fusion').ViewObject.ShapeColor)
App.getDocument('Combustor').getObject('Fusion').ViewObject.DisplayMode=getattr(App.getDocument('Combustor').getObject('Cut002').getLinkedObject(True).ViewObject,'DisplayMode',App.getDocument('Combustor').getObject('Fusion').ViewObject.DisplayMode)
App.ActiveDocument.recompute()

###########
############
#############
###########
############
############

# Now generate flames again
App.ActiveDocument.addObject("Part::Cylinder","Cylinder")
App.ActiveDocument.ActiveObject.Label = "Cylinder010"
App.ActiveDocument.recompute()
FreeCAD.getDocument('Combustor').getObject('Cylinder010').Radius = R_flame.str
FreeCAD.getDocument('Combustor').getObject('Cylinder010').Height = L_flame.str
Gui.SendMsgToActiveView("ViewFit")

FreeCAD.getDocument('Combustor').getObject('Cylinder010').Placement = App.Placement(App.Vector(R_mid_plane.mm,0,0),App.Rotation(App.Vector(0,0,1),0))

_obj_ = Draft.make_polar_array(App.ActiveDocument.Cylinder010, number=N_s, angle=360.0, center=FreeCAD.Vector(0.0, 0.0, 0.0), use_link=True)
_obj_.Fuse = False
Draft.autogroup(_obj_)
App.ActiveDocument.recompute()

##
#
#
#
#
#
#
# Now export step files

flame_name = path+"flame.step"
fusion_name = path+"fusion.step"

__objs__=[]
__objs__.append(FreeCAD.getDocument("Combustor").getObject("Fusion"))
import ImportGui
ImportGui.export(__objs__,fusion_name)

del __objs__

__objs__=[]
__objs__.append(FreeCAD.getDocument("Combustor").getObject("Array004"))
ImportGui.export(__objs__,flame_name)

del __objs__
