"""
cbunge
Washington State University
email: carl.bunge@wsu.edu

Adapted from @author: Luka Denies from TU Delft.

Changelog:
cbunge - 11/2017 - Integration of CoolProp
cbunge - 06/2018 - Update to OpenFOAM-5.x (Mass-based thermodynamics (for example: cpMcv to CpMCv))
cbunge - 07/2018 - Implemented TTable lookup format into heTabularThermo class. No need to directly alter core thermo files with an inverse interpolate function which original architecture required. Be sure to check the TTable reference property for your application. For sensibleInternalEnergy  use the E (internal energy) implentation. For sensibleEnthalpy use H. See note above TTable write script below.

"""

import CoolProp.CoolProp as CP
import numpy as np
import matplotlib.pyplot as plt

#Fluid for thermodynamic properties (rho, Cp, CpMcv, H)
fluid_thermo ='parahydrogen'

#Fluid for transport model (thermal conductivity and viscosity)
fluid_transport = 'hydrogen'

#Universal gas constant, R [J/(mol K)]
#R = CP.PropsSI("GAS_CONSTANT",fluid_thermo)

#****************************************************************************************

#Temperature limits
T0 = 19.3 # - 19.3K is 0.03K above saturated conditions (@75kPa)
TMax = 75 #K

#Pressure limits
p0 = 7.5e4 #Pa
pMax = 10e5 #Pa

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
E = []

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
    E.append([0])
    rho[i][0] = rhoCur = CP.PropsSI('D','T',T,'P',p,fluid_thermo)
    CpCur = CP.PropsSI('C','D',rhoCur,'T',T,fluid_thermo) 
    Cp[i][0] = CpCur
    mu[i][0] = CP.PropsSI('V','D',rhoCur,'T',T,fluid_transport)
    kappa[i][0] = CP.PropsSI('L','D',rhoCur,'T',T,fluid_transport) 
    CpMCv[i][0] = CpCur-CP.PropsSI('O','D',rhoCur,'T',T,fluid_thermo) 
    H[i][0] =  CP.PropsSI('H','D',rhoCur,'T',T,fluid_thermo)
    E[i][0] =  CP.PropsSI('U','D',rhoCur,'T',T,fluid_thermo)
    TRange.append(T)
    while T<TMax:
        j += 1
        dT = 0.1 # T [K] **************************************************************
        T += dT
        rhoCur = CP.PropsSI('D','T',T,'P',p,fluid_thermo)
        rho[i].append(rhoCur)
        CpCur = CP.PropsSI('C','D',rhoCur,'T',T,fluid_thermo)
        Cp[i].append(CpCur)
        mu[i].append(CP.PropsSI('V','D',rhoCur,'T',T,fluid_transport))
        kappa[i].append(CP.PropsSI('L','D',rhoCur,'T',T,fluid_transport))
        CpMCv[i].append((CpCur-CP.PropsSI('O','D',rhoCur,'T',T,fluid_thermo)))
        H[i].append(CP.PropsSI('H','D',rhoCur,'T',T,fluid_thermo))
        E[i].append(CP.PropsSI('U','D',rhoCur,'T',T,fluid_thermo)) 
        TRange.append(T)
    i += 1
    ps.append([p]*len(TRange))    
    rhoPseudoCrit = CP.PropsSI('D','T',Tcrit,'P',p,fluid_thermo)
    dp = 1000 # P [Pa] ****************************************************************
    p += dp
    print p
    Ts.append(TRange)
print "Calculations done, now writing"

muFile = open("mu","w")
muFile.write("( \n")

for i,p in enumerate(pRange):
    muFile.write("(" + str(p) + "\n(\n")
    sList = ["\t(" + str(Ts[i][j]) + " " + str(mu[i][j]) + ")\n" for j in range(len(Ts[i]))]
    muFile.write(" ".join(sList))
    muFile.write(") ) \n")    
muFile.write(");")
muFile.close()

rhoFile = open("rho","w")
rhoFile.write("( \n")

for i,p in enumerate(pRange):
    rhoFile.write("(" + str(p) + "\n(\n")
    sList = ["\t(" + str(Ts[i][j]) + " " + str(rho[i][j]) + ")\n" for j in range(len(Ts[i]))]
    rhoFile.write(" ".join(sList))
    rhoFile.write(") ) \n")    
rhoFile.write(");")
rhoFile.close()

CpFile = open("Cp","w")
CpFile.write("( \n")

for i,p in enumerate(pRange):
    CpFile.write("(" + str(p) + "\n(\n")
    sList = ["\t(" + str(Ts[i][j]) + " " + str(Cp[i][j]) + ")\n" for j in range(len(Ts[i]))]
    CpFile.write(" ".join(sList))
    CpFile.write(") ) \n")    
CpFile.write(");")
CpFile.close()

kappaFile = open("kappa","w")
kappaFile.write("( \n")

for i,p in enumerate(pRange):
    kappaFile.write("(" + str(p) + "\n(\n")
    sList = ["\t(" + str(Ts[i][j]) + " " + str(kappa[i][j]) + ")\n" for j in range(len(Ts[i]))]
    kappaFile.write(" ".join(sList))
    kappaFile.write(") ) \n")    
kappaFile.write(");")
kappaFile.close()

CpMCvFile = open("CpMCv","w")
CpMCvFile.write("( \n")

for i,p in enumerate(pRange):
    CpMCvFile.write("(" + str(p) + "\n(\n")
    sList = ["\t(" + str(Ts[i][j]) + " " + str(CpMCv[i][j]) + ")\n" for j in range(len(Ts[i]))]
    CpMCvFile.write(" ".join(sList))
    CpMCvFile.write(") ) \n")    
CpMCvFile.write(");")
CpMCvFile.close()

HFile = open("H","w")
HFile.write("( \n")

for i,p in enumerate(pRange):
    HFile.write("(" + str(p) + "\n(\n")
    sList = ["\t(" + str(Ts[i][j]) + " " + str(H[i][j]) + ")\n" for j in range(len(Ts[i]))]
    HFile.write(" ".join(sList))
    HFile.write(") ) \n")    
HFile.write(");")
HFile.close()

EFile = open("E","w")
EFile.write("( \n")

for i,p in enumerate(pRange):
    EFile.write("(" + str(p) + "\n(\n")
    sList = ["\t(" + str(Ts[i][j]) + " " + str(E[i][j]) + ")\n" for j in range(len(Ts[i]))]
    EFile.write(" ".join(sList))
    EFile.write(") ) \n")    
EFile.write(");")
EFile.close()

## TTable - Currrently InternalEnergy (E) -based (change E to H for Enthalpy based temperature lookup ##
TFile = open("TTable","w")
TFile.write("( \n")

for i,p in enumerate(pRange):
    TFile.write("(" + str(p) + "\n(\n")
    sList = ["\t(" + str(E[i][j]) + " " + str(Ts[i][j]) + ")\n" for j in range(len(Ts[i]))]
    TFile.write(" ".join(sList))
    TFile.write(") ) \n")    
TFile.write(");")
TFile.close()


#Previous dT method to save computational time (was not able to produce desired outcome):
#dT = drho/CP.PropsSI('d(D)/d(P)|T','D',rhoCur,'T',T,fluid_thermo)*CP.PropsSI('d(P)/d(T)|D','D',rhoCur,'T',T,fluid_thermo)

#Previous dP method to save computational time (was not able to produce desired outcome):
#drho/CP.PropsSI('d(D)/d(P)|T','D',rhoPseudoCrit,'T',Tcrit,fluid_thermo)
