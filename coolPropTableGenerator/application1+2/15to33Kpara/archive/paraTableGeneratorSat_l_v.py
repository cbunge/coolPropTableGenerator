"""
Carl Bunge
Washington State University
June 2018

Adapted from @author: Luka Denies from TU Delft.

Changelog:
11/2017 - Integration of CoolProp
06/2018 - Update to OpenFOAM-5.x (Mass-based thermodynamics (for example: cpMcv to CpMCv))
03/2019 - Update to include parahydrogen properties from Refprop

"""

import CoolProp.CoolProp as CP
import numpy as np
import matplotlib.pyplot as plt

#Fluid for thermodynamic properties (rho, Cp, CpMcv, H, S, c, E, thermal conductivity)
CP.set_reference_state('parahydrogen','NBP')
fluid_thermo ='parahydrogen'

#Fluid for transport model (viscosity)
#CP.set_reference_state('hydrogen','NBP')
fluid_transport = 'hydrogen'

#****************************************************************************************

#Temperature limits (set within subcritical region for saturation tables)
T0 = 15 #Temperature start (K)
TMax = 32 #Temperature end (K)

#Pressure limits
p0 = 0.1e5 #Pa
pMax = 5.5e5 #Pa

#****************************************************************************************

Tcrit = CP.PropsSI("Tcrit",fluid_thermo)

Ts = []
ps = []
pRange = []

rho = []
mu = []
mu_l = []
mu_v = []
kappa = []
kappa_l = []
kappa_v = []
Cp = []
Cp_l = []
Cp_v = []
H = []
H_l = []
H_v = []
CpMCv = []
CpMCv_l = []
CpMCv_v = []
E = []
E_l = []
E_v = []
S = []
S_l = []
S_v = []
c = []
c_l = []
c_v = []
pSat = []

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
    Cp_l.append([0])
    Cp_v.append([0])
    mu.append([0])
    mu_l.append([0])
    mu_v.append([0])
    kappa.append([0])
    kappa_l.append([0])
    kappa_v.append([0])
    CpMCv.append([0])
    CpMCv_l.append([0])
    CpMCv_v.append([0])
    H.append([0])
    H_l.append([0])
    H_v.append([0])
    E.append([0])
    E_l.append([0])
    E_v.append([0])
    S.append([0])
    S_l.append([0])
    S_v.append([0])
    c.append([0])
    c_l.append([0])
    c_v.append([0])
    pSat.append([0])
    rho[i][0] = rhoCur = CP.PropsSI('D','T',T,'P',p,fluid_thermo)
    CpCur = CP.PropsSI('C','D',rhoCur,'T',T,fluid_thermo) 
    Cp[i][0] = CpCur
    Cp_l[i][0] = CP.PropsSI('C','T',T,'Q',0,fluid_thermo) 
    Cp_v[i][0] = CP.PropsSI('C','T',T,'Q',1,fluid_thermo)
    mu_l[i][0] = CP.PropsSI('V','T',T,'Q',0,fluid_transport) 
    mu_v[i][0] = CP.PropsSI('V','T',T,'Q',1,fluid_transport)
    mu[i][0] = CP.PropsSI('V','D',rhoCur,'T',T,fluid_transport)
    kappa_l[i][0] = CP.PropsSI('L','T',T,'Q',0,'REFPROP::parahydrogen') 
    kappa_v[i][0] = CP.PropsSI('L','T',T,'Q',1,'REFPROP::parahydrogen')
    kappa[i][0] = CP.PropsSI('L','D',rhoCur,'T',T,'REFPROP::parahydrogen') 
    CpMCv_l[i][0] = CP.PropsSI('O','T',T,'Q',0,fluid_thermo) 
    CpMCv_v[i][0] = CP.PropsSI('O','T',T,'Q',1,fluid_thermo)
    CpMCv[i][0] = CpCur-CP.PropsSI('O','D',rhoCur,'T',T,fluid_thermo) 
    H_l[i][0] = CP.PropsSI('H','T',T,'Q',0,fluid_thermo) 
    H_v[i][0] = CP.PropsSI('H','T',T,'Q',1,fluid_thermo)
    H[i][0] =  CP.PropsSI('H','D',rhoCur,'T',T,fluid_thermo)
    E_l[i][0] = CP.PropsSI('U','T',T,'Q',0,fluid_thermo) 
    E_v[i][0] = CP.PropsSI('U','T',T,'Q',1,fluid_thermo)
    E[i][0] =  CP.PropsSI('U','D',rhoCur,'T',T,fluid_thermo) 
    S_l[i][0] = CP.PropsSI('S','T',T,'Q',0,fluid_thermo) 
    S_v[i][0] = CP.PropsSI('S','T',T,'Q',1,fluid_thermo)
    S[i][0] =  CP.PropsSI('S','D',rhoCur,'T',T,fluid_thermo)
    c_l[i][0] = CP.PropsSI('A','T',T,'Q',0,fluid_thermo) 
    c_v[i][0] = CP.PropsSI('A','T',T,'Q',1,fluid_thermo) 
    c[i][0] =  CP.PropsSI('A','D',rhoCur,'T',T,fluid_thermo)
    pSat[i][0] = CP.PropsSI('P','T',T,'Q',0,fluid_thermo)
    TRange.append(T)
    while T<TMax:
        j += 1
        dT = 1 # Tstep [K] **************************************************************
        T += dT
        rhoCur = CP.PropsSI('D','T',T,'P',p,fluid_thermo)
        rho[i].append(rhoCur)
        CpCur = CP.PropsSI('C','D',rhoCur,'T',T,fluid_thermo) 
        CpCur_l = CP.PropsSI('C','T',T,'Q',0,fluid_thermo)
        CpCur_v = CP.PropsSI('C','T',T,'Q',1,fluid_thermo) 
        Cp_l[i].append(CP.PropsSI('C','T',T,'Q',0,fluid_thermo))
        Cp_v[i].append(CP.PropsSI('C','T',T,'Q',1,fluid_thermo))
        Cp[i].append(CpCur)
        mu_l[i].append(CP.PropsSI('V','T',T,'Q',0,fluid_transport))
        mu_v[i].append(CP.PropsSI('V','T',T,'Q',1,fluid_transport))
        mu[i].append(CP.PropsSI('V','D',rhoCur,'T',T,fluid_transport))
        kappa_l[i].append(CP.PropsSI('L','T',T,'Q',0,'REFPROP::parahydrogen'))
        kappa_v[i].append(CP.PropsSI('L','T',T,'Q',1,'REFPROP::parahydrogen'))
        kappa[i].append(CP.PropsSI('L','D',rhoCur,'T',T,'REFPROP::parahydrogen'))
        CpMCv_l[i].append((CP.PropsSI('C','T',T,'Q',0,fluid_thermo))-(CP.PropsSI('O','T',T,'Q',0,fluid_thermo)))
        CpMCv_v[i].append((CP.PropsSI('C','T',T,'Q',1,fluid_thermo))-(CP.PropsSI('O','T',T,'Q',1,fluid_thermo)))
        CpMCv[i].append((CpCur-CP.PropsSI('O','D',rhoCur,'T',T,fluid_thermo)))
        H_l[i].append(CP.PropsSI('H','T',T,'Q',0,fluid_thermo))
        H_v[i].append(CP.PropsSI('H','T',T,'Q',1,fluid_thermo))
        H[i].append(CP.PropsSI('H','D',rhoCur,'T',T,fluid_thermo))
        E_l[i].append(CP.PropsSI('U','T',T,'Q',0,fluid_thermo))
        E_v[i].append(CP.PropsSI('U','T',T,'Q',1,fluid_thermo))
        E[i].append(CP.PropsSI('U','D',rhoCur,'T',T,fluid_thermo))
        S_l[i].append(CP.PropsSI('S','T',T,'Q',0,fluid_thermo))
        S_v[i].append(CP.PropsSI('S','T',T,'Q',1,fluid_thermo))
        S[i].append(CP.PropsSI('S','D',rhoCur,'T',T,fluid_thermo))
        c_l[i].append(CP.PropsSI('A','T',T,'Q',0,fluid_thermo))
        c_v[i].append(CP.PropsSI('A','T',T,'Q',1,fluid_thermo))
        c[i].append(CP.PropsSI('A','D',rhoCur,'T',T,fluid_thermo))
        pSat[i].append(CP.PropsSI('P','T',T,'Q',0,fluid_thermo))
        TRange.append(T)
    i += 1
    ps.append([p]*len(TRange))    
    rhoPseudoCrit = CP.PropsSI('D','T',Tcrit,'P',p,fluid_thermo)
    dp = 0.5e5 # Pstep [Pa] ****************************************************************
    p += dp
    print p
    Ts.append(TRange)
print "Calculations done, now writing"

pSatFile = open("pSat","w")

for i,p in enumerate(pRange):
    sList = ["\t" + str(pSat[i][j]) + " " + str(Ts[i][j]) + "\n" for j in range(len(Ts[i]))]
    pSatFile.write("".join(sList))    
pSatFile.write("")
pSatFile.close()

mu_lFile = open("mu_l","w")

for i,p in enumerate(pRange):
    sList = ["\t" + str(mu_l[i][j]) + " " + str(Ts[i][j]) + " " +  str(p) + "\n" for j in range(len(Ts[i]))]
    mu_lFile.write("".join(sList))    
mu_lFile.write("")
mu_lFile.close()

mu_vFile = open("mu_v","w")

for i,p in enumerate(pRange):
    sList = ["\t" + str(mu_v[i][j]) + " " + str(Ts[i][j]) + " " +  str(p) + "\n" for j in range(len(Ts[i]))]
    mu_vFile.write("".join(sList))    
mu_vFile.write("")
mu_vFile.close()

muFile = open("mu","w")

for i,p in enumerate(pRange):
    sList = ["\t" + str(mu[i][j]) + " " + str(Ts[i][j]) + " " +  str(p) + "\n" for j in range(len(Ts[i]))]
    muFile.write("".join(sList))    
muFile.write("")
muFile.close()

rhoFile = open("rho","w")
rhoFile.write("\n")

for i,p in enumerate(pRange):
    rhoFile.write("")
    sList = ["\t" + str(rho[i][j]) + " " + str(Ts[i][j]) + " " +  str(p) + "\n" for j in range(len(Ts[i]))]
    rhoFile.write("".join(sList))
rhoFile.write("")
rhoFile.close()

Cp_lFile = open("Cp_l","w")
Cp_lFile.write("\n")

for i,p in enumerate(pRange):
    Cp_lFile.write("")
    sList = ["\t" + str(Cp_l[i][j]) + " " + str(Ts[i][j]) + " " +  str(p) + "\n" for j in range(len(Ts[i]))]
    Cp_lFile.write("".join(sList))
Cp_lFile.write("")
Cp_lFile.close()

Cp_vFile = open("Cp_v","w")
Cp_vFile.write("\n")

for i,p in enumerate(pRange):
    Cp_vFile.write("")
    sList = ["\t" + str(Cp_v[i][j]) + " " + str(Ts[i][j]) + " " +  str(p) + "\n" for j in range(len(Ts[i]))]
    Cp_vFile.write("".join(sList))
Cp_vFile.write("")
Cp_vFile.close()

CpFile = open("Cp","w")
CpFile.write("\n")

for i,p in enumerate(pRange):
    CpFile.write("")
    sList = ["\t" + str(Cp[i][j]) + " " + str(Ts[i][j]) + " " +  str(p) + "\n" for j in range(len(Ts[i]))]
    CpFile.write("".join(sList))
CpFile.write("")
CpFile.close()

kappa_lFile = open("kappa_l","w")
kappa_lFile.write("\n")

for i,p in enumerate(pRange):
    kappa_lFile.write("")
    sList = ["\t" + str(kappa_l[i][j]) + " " + str(Ts[i][j]) + " " +  str(p) + "\n" for j in range(len(Ts[i]))]
    kappa_lFile.write("".join(sList))
kappa_lFile.write("")
kappa_lFile.close()

kappa_vFile = open("kappa_v","w")
kappa_vFile.write("\n")

for i,p in enumerate(pRange):
    kappa_vFile.write("")
    sList = ["\t" + str(kappa_v[i][j]) + " " + str(Ts[i][j]) + " " +  str(p) + "\n" for j in range(len(Ts[i]))]
    kappa_vFile.write("".join(sList))
kappa_vFile.write("")
kappa_vFile.close()

kappaFile = open("kappa","w")
kappaFile.write("\n")

for i,p in enumerate(pRange):
    kappaFile.write("")
    sList = ["\t" + str(kappa[i][j]) + " " + str(Ts[i][j]) + " " +  str(p) + "\n" for j in range(len(Ts[i]))]
    kappaFile.write("".join(sList))
kappaFile.write("")
kappaFile.close()

CpMCv_lFile = open("CpMCv_l","w")
CpMCv_lFile.write("\n")

for i,p in enumerate(pRange):
    CpMCv_lFile.write("")
    sList = ["\t" + str(CpMCv_l[i][j]) + " " + str(Ts[i][j]) + " " +  str(p) + "\n" for j in range(len(Ts[i]))]
    CpMCv_lFile.write("".join(sList))
CpMCv_lFile.write("")
CpMCv_lFile.close()

CpMCv_vFile = open("CpMCv_v","w")
CpMCv_vFile.write("\n")

for i,p in enumerate(pRange):
    CpMCv_vFile.write("")
    sList = ["\t" + str(CpMCv_v[i][j]) + " " + str(Ts[i][j]) + " " +  str(p) + "\n" for j in range(len(Ts[i]))]
    CpMCv_vFile.write("".join(sList))
CpMCv_vFile.write("")
CpMCv_vFile.close()

CpMCvFile = open("CpMCv","w")
CpMCvFile.write("\n")

for i,p in enumerate(pRange):
    CpMCvFile.write("")
    sList = ["\t" + str(CpMCv[i][j]) + " " + str(Ts[i][j]) + " " +  str(p) + "\n" for j in range(len(Ts[i]))]
    CpMCvFile.write("".join(sList))
CpMCvFile.write("")
CpMCvFile.close()

H_lFile = open("H_l","w")
H_lFile.write("\n")

for i,p in enumerate(pRange):
    H_lFile.write("")
    sList = ["\t" + str(H_l[i][j]) + " " + str(Ts[i][j]) + " " +  str(p) + "\n" for j in range(len(Ts[i]))]
    H_lFile.write("".join(sList))
H_lFile.write("")
H_lFile.close()

H_vFile = open("H_v","w")
H_vFile.write("\n")

for i,p in enumerate(pRange):
    H_vFile.write("")
    sList = ["\t" + str(H_v[i][j]) + " " + str(Ts[i][j]) + " " +  str(p) + "\n" for j in range(len(Ts[i]))]
    H_vFile.write("".join(sList))
H_vFile.write("")
H_vFile.close()

HFile = open("H","w")
HFile.write("\n")

for i,p in enumerate(pRange):
    HFile.write("")
    sList = ["\t" + str(H[i][j]) + " " + str(Ts[i][j]) + " " +  str(p) + "\n" for j in range(len(Ts[i]))]
    HFile.write("".join(sList))
HFile.write("")
HFile.close()

E_lFile = open("E_l","w")
E_lFile.write("\n")

for i,p in enumerate(pRange):
    E_lFile.write("")
    sList = ["\t" + str(E_l[i][j]) + " " + str(Ts[i][j]) + " " +  str(p) + "\n" for j in range(len(Ts[i]))]
    E_lFile.write("".join(sList))
E_lFile.write("")
E_lFile.close()

E_vFile = open("E_v","w")
E_vFile.write("\n")

for i,p in enumerate(pRange):
    E_vFile.write("")
    sList = ["\t" + str(E_v[i][j]) + " " + str(Ts[i][j]) + " " +  str(p) + "\n" for j in range(len(Ts[i]))]
    E_vFile.write("".join(sList))
E_vFile.write("")
E_vFile.close()

EFile = open("E","w")
EFile.write("\n")

for i,p in enumerate(pRange):
    EFile.write("")
    sList = ["\t" + str(E[i][j]) + " " + str(Ts[i][j]) + " " +  str(p) + "\n" for j in range(len(Ts[i]))]
    EFile.write("".join(sList))
EFile.write("")
EFile.close()

S_lFile = open("S_l","w")
S_lFile.write("\n")

for i,p in enumerate(pRange):
    S_lFile.write("")
    sList = ["\t" + str(S_l[i][j]) + " " + str(Ts[i][j]) + " " +  str(p) + "\n" for j in range(len(Ts[i]))]
    S_lFile.write("".join(sList))
S_lFile.write("")
S_lFile.close()

S_vFile = open("S_v","w")
S_vFile.write("\n")

for i,p in enumerate(pRange):
    S_vFile.write("")
    sList = ["\t" + str(S_v[i][j]) + " " + str(Ts[i][j]) + " " +  str(p) + "\n" for j in range(len(Ts[i]))]
    S_vFile.write("".join(sList))
S_vFile.write("")
S_vFile.close()

SFile = open("S","w")
SFile.write("\n")

for i,p in enumerate(pRange):
    SFile.write("")
    sList = ["\t" + str(S[i][j]) + " " + str(Ts[i][j]) + " " +  str(p) + "\n" for j in range(len(Ts[i]))]
    SFile.write("".join(sList))
SFile.write("")
SFile.close()

c_lFile = open("c_l","w")
c_lFile.write("\n")

for i,p in enumerate(pRange):
    c_lFile.write("")
    sList = ["\t" + str(c_l[i][j]) + " " + str(Ts[i][j]) + " " +  str(p) + "\n" for j in range(len(Ts[i]))]
    c_lFile.write("".join(sList))
c_lFile.write("")
c_lFile.close()

c_vFile = open("c_v","w")
c_vFile.write("\n")

for i,p in enumerate(pRange):
    c_vFile.write("")
    sList = ["\t" + str(c_v[i][j]) + " " + str(Ts[i][j]) + " " +  str(p) + "\n" for j in range(len(Ts[i]))]
    c_vFile.write("".join(sList))
c_vFile.write("")
c_vFile.close()

cFile = open("c","w")
cFile.write("\n")

for i,p in enumerate(pRange):
    cFile.write("")
    sList = ["\t" + str(c[i][j]) + " " + str(Ts[i][j]) + " " +  str(p) + "\n" for j in range(len(Ts[i]))]
    cFile.write("".join(sList))
cFile.write("")
cFile.close()

#Previous dT method to save computational time:
#dT = drho/CP.PropsSI('d(D)/d(P)|T','D',rhoCur,'T',T,fluid_thermo)*CP.PropsSI('d(P)/d(T)|D','D',rhoCur,'T',T,fluid_thermo)

#Previous dP method to save computational time:
#drho/CP.PropsSI('d(D)/d(P)|T','D',rhoPseudoCrit,'T',Tcrit,fluid_thermo)
