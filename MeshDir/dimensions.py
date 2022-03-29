

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
