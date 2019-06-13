"""
Carl Bunge
Washington State University
June 2018

Adapted from @author: Luka Denies from TU Delft.

Changelog:
11/2017 - Integration of CoolProp
06/2018 - Update to OpenFOAM-5.x (Mass-based thermodynamics (for example: cpMcv to CpMCv))

"""

import CoolProp.CoolProp as CP
import numpy as np
import matplotlib.pyplot as plt

#Fluid for thermodynamic properties (rho, Cp, CpMcv, H)
fluid_thermo ='orthohydrogen'

#Fluid for transport model (thermal conductivity and viscosity)
fluid_transport = 'hydrogen'

#Universal gas constant, R [J/(mol K)]
#R = CP.PropsSI("GAS_CONSTANT",fluid_thermo)

#****************************************************************************************

#Temperature limits
T0 = 20 # - 19.3K is 0.03K above saturated conditions
TMax = 90 #

#Pressure limits
p0 = 0.5e5 #Pa
pMax = 5.5e5 #Pa

#****************************************************************************************

Tcrit = CP.PropsSI("Tcrit",fluid_thermo)

Ts = []
ps = []
pRange = []

rho = []
mu = []
kappa = []
Cp = []
H = []
CpMCv = []

i = 0
j = 0

p = p0
T = T0

#Build (p, T) tables
while p<pMax:
    pRange.append(p)
    TRange = []
    T = T0
    rho.append([0])
    Cp.append([0])
    mu.append([0])
    kappa.append([0])
    CpMCv.append([0])
    H.append([0])
    rho[i][0] = rhoCur = CP.PropsSI('D','T',T,'P',p,fluid_thermo)
    CpCur = CP.PropsSI('C','D',rhoCur,'T',T,fluid_thermo) 
    Cp[i][0] = CpCur
    mu[i][0] = CP.PropsSI('V','D',rhoCur,'T',T,fluid_transport)
    kappa[i][0] = CP.PropsSI('L','D',rhoCur,'T',T,fluid_transport) 
    CpMCv[i][0] = CpCur-CP.PropsSI('O','D',rhoCur,'T',T,fluid_thermo) 
    H[i][0] =  CP.PropsSI('H','D',rhoCur,'T',T,fluid_thermo)
    TRange.append(T)
    while T<TMax:
        j += 1
        dT = 10 # Tstep [K] **************************************************************
        T += dT
        rhoCur = CP.PropsSI('D','T',T,'P',p,fluid_thermo)
        rho[i].append(rhoCur)
        CpCur = CP.PropsSI('C','D',rhoCur,'T',T,fluid_thermo)
        Cp[i].append(CpCur)
        mu[i].append(CP.PropsSI('V','D',rhoCur,'T',T,fluid_transport))
        kappa[i].append(CP.PropsSI('L','D',rhoCur,'T',T,fluid_transport))
        CpMCv[i].append((CpCur-CP.PropsSI('O','D',rhoCur,'T',T,fluid_thermo)))
        H[i].append(CP.PropsSI('H','D',rhoCur,'T',T,fluid_thermo))
        TRange.append(T)
    i += 1
    ps.append([p]*len(TRange))    
    rhoPseudoCrit = CP.PropsSI('D','T',Tcrit,'P',p,fluid_thermo)
    dp = 0.5e5 # Pstep [Pa] ****************************************************************
    p += dp
    print p
    Ts.append(TRange)
print "Calculations done, now writing"

muFile = open("mu","w")
muFile.write("\n")

for i,p in enumerate(pRange):
    muFile.write("")
    sList = ["\t" + str(mu[i][j]) + " " + str(Ts[i][j]) + " " +  str(p) + "\n" for j in range(len(Ts[i]))]
    muFile.write(" ".join(sList))
    muFile.write("")    
muFile.write("")
muFile.close()

rhoFile = open("rho","w")
rhoFile.write("\n")

for i,p in enumerate(pRange):
    rhoFile.write("")
    sList = ["\t" + str(rho[i][j]) + " " + str(Ts[i][j]) + " " +  str(p) + "\n" for j in range(len(Ts[i]))]
    rhoFile.write(" ".join(sList))
    rhoFile.write(" \n")    
rhoFile.write("")
rhoFile.close()

CpFile = open("Cp","w")
CpFile.write("\n")

for i,p in enumerate(pRange):
    CpFile.write("")
    sList = ["\t" + str(Cp[i][j]) + " " + str(Ts[i][j]) + " " +  str(p) + "\n" for j in range(len(Ts[i]))]
    CpFile.write(" ".join(sList))
    CpFile.write(" \n")    
CpFile.write("")
CpFile.close()

kappaFile = open("kappa","w")
kappaFile.write("\n")

for i,p in enumerate(pRange):
    kappaFile.write("")
    sList = ["\t" + str(kappa[i][j]) + " " + str(Ts[i][j]) + " " +  str(p) + "\n" for j in range(len(Ts[i]))]
    kappaFile.write(" ".join(sList))
    kappaFile.write(" \n")    
kappaFile.write("")
kappaFile.close()

CpMCvFile = open("CpMCv","w")
CpMCvFile.write("\n")

for i,p in enumerate(pRange):
    CpMCvFile.write("")
    sList = ["\t" + str(CpMCv[i][j]) + " " + str(Ts[i][j]) + " " +  str(p) + "\n" for j in range(len(Ts[i]))]
    CpMCvFile.write(" ".join(sList))
    CpMCvFile.write(" \n")    
CpMCvFile.write("")
CpMCvFile.close()

HFile = open("H","w")
HFile.write("\n")

for i,p in enumerate(pRange):
    HFile.write("")
    sList = ["\t" + str(H[i][j]) + " " + str(Ts[i][j]) + " " +  str(p) + "\n" for j in range(len(Ts[i]))]
    HFile.write(" ".join(sList))
    HFile.write(" \n")    
HFile.write("")
HFile.close()

#Previous dT method to save computational time:
#dT = drho/CP.PropsSI('d(D)/d(P)|T','D',rhoCur,'T',T,fluid_thermo)*CP.PropsSI('d(P)/d(T)|D','D',rhoCur,'T',T,fluid_thermo)

#Previous dP method to save computational time:
#drho/CP.PropsSI('d(D)/d(P)|T','D',rhoPseudoCrit,'T',Tcrit,fluid_thermo)
