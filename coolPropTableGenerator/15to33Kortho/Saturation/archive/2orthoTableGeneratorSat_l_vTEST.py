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

#Custom refernce state for orthohydrogen to capture enthalpy offset of orthohydrogen
CP.set_reference_state('orthohydrogen',20.3800689304,35150.6373702,1417.12332,0.036828) 
fluid_thermo ='orthohydrogen'

#Fluid for transport model (viscosity)
fluid_transport = 'hydrogen'

#****************************************************************************************

#Temperature limits (set within subcritical region for saturation tables)
T0 = 15 #Temperature start (K)
TMax = 32 #Temperature end (K)

#Pressure limits
p0 = CP.PropsSI('P','T',T0,'Q',0,fluid_thermo) #Pa
pMax = CP.PropsSI('P','T',TMax,'Q',0,fluid_thermo) #Pa

#****************************************************************************************

Tcrit = CP.PropsSI("Tcrit",fluid_thermo)

Ts = []
pSats = []
TRange = []

rho_l = []
rho_v = []
mu_l = []
mu_v = []
kappa_l = []
kappa_v = []
Cp_l = []
Cp_v = []
H_l = []
H_v = []
CpMCv_l = []
CpMCv_v = []
E_l = []
E_v = []
S_l = []
S_v = []
c_l = []
c_v = []
pSat = []

i = 0
j = 0

p = p0
T = T0-1

#Build (p, T) tables
while T<TMax:      
    TRange.append(T)
    pSatRange = []
    
    rho_l.append([0])
    rho_v.append([0])
    Cp_l.append([0])
    Cp_v.append([0])
    mu_l.append([0])
    mu_v.append([0])
    kappa_l.append([0])
    kappa_v.append([0])
    CpMCv_l.append([0])
    CpMCv_v.append([0])
    H_l.append([0])
    H_v.append([0])
    E_l.append([0])
    E_v.append([0])
    S_l.append([0])
    S_v.append([0])
    c_l.append([0])
    c_v.append([0])
    pSat.append([0])

    pSat[i].append(CP.PropsSI('P','T',T,'Q',0,fluid_thermo)) 
    rho_l[i].append(CP.PropsSI('D','T',T,'Q',0,fluid_thermo))
    rho_v[i].append(CP.PropsSI('D','T',T,'Q',1,fluid_thermo))
    Cp_l[i].append(CP.PropsSI('C','T',T,'Q',0,fluid_thermo))
    Cp_v[i].append(CP.PropsSI('C','T',T,'Q',1,fluid_thermo))
    mu_l[i].append(CP.PropsSI('V','T',T,'Q',0,fluid_transport))
    mu_v[i].append(CP.PropsSI('V','T',T,'Q',1,fluid_transport))
    kappa_l[i].append((CP.PropsSI('L','T',T,'Q',0,fluid_transport)-(0.25*CP.PropsSI('L','T',T,'Q',0,'REFPROP::parahydrogen')))*(1.3333333))
    kappa_v[i].append((CP.PropsSI('L','T',T,'Q',1,fluid_transport)-(0.25*CP.PropsSI('L','T',T,'Q',1,'REFPROP::parahydrogen')))*(1.3333333))
    CpMCv_l[i].append((CP.PropsSI('C','T',T,'Q',0,fluid_thermo))-(CP.PropsSI('O','T',T,'Q',0,fluid_thermo)))
    CpMCv_v[i].append((CP.PropsSI('C','T',T,'Q',1,fluid_thermo))-(CP.PropsSI('O','T',T,'Q',1,fluid_thermo)))
    H_l[i].append(CP.PropsSI('H','T',T,'Q',0,fluid_thermo))
    H_v[i].append(CP.PropsSI('H','T',T,'Q',1,fluid_thermo))
    E_l[i].append(CP.PropsSI('U','T',T,'Q',0,fluid_thermo))
    E_v[i].append(CP.PropsSI('U','T',T,'Q',1,fluid_thermo))
    S_l[i].append(CP.PropsSI('S','T',T,'Q',0,fluid_thermo))
    S_v[i].append(CP.PropsSI('S','T',T,'Q',1,fluid_thermo))
    c_l[i].append(CP.PropsSI('A','T',T,'Q',0,fluid_thermo))
    c_v[i].append(CP.PropsSI('A','T',T,'Q',1,fluid_thermo))
    pSatRange.append(pSat)
    
    i += 1
    dT = 1 # Tstep [K] **************************************************************
    T += dT 
    
    Ts.append([T]*len(pSatRange))         
    pSats.append(pSatRange)
    print T
print "Calculations done, now writing"

mu_lFile = open("mu_l","w")

for i,T in enumerate(TRange):
    raw_mu_l = str(mu_l[i])[1:-1]
    mu_l1 = raw_mu_l.replace("0","")
    mu_l2 = mu_l1.replace(" ", "")
    mu_l3 = mu_l2.replace(",","")
    raw_pSat = str(pSat[i])[1:-1]
    pSat_l1 = raw_pSat.replace("0","")
    pSat_l2 = pSat_l1.replace(" ", "")
    pSat_l3 = pSat_l2.replace(",","")
    
    sList = ["\t" + mu_l3 + " " + str(Ts[i][j]) + " " +  pSat_l3 + "\n" for j in range(len(Ts[i]))]
    mu_lFile.write("".join(sList))    
mu_lFile.write("")
mu_lFile.close()

mu_vFile = open("mu_v","w")

for i,T in enumerate(TRange):
    raw_mu_v = str(mu_v[i])[1:-1]
    mu_v1 = raw_mu_v.replace("0","")
    mu_v2 = mu_v1.replace(" ", "")
    mu_v3 = mu_v2.replace(",","")
    raw_pSat = str(pSat[i])[1:-1]
    pSat_l1 = raw_pSat.replace("0","")
    pSat_l2 = pSat_l1.replace(" ", "")
    pSat_l3 = pSat_l2.replace(",","")
    
    sList = ["\t" + mu_v3 + " " + str(Ts[i][j]) + " " + pSat_l3 + "\n" for j in range(len(Ts[i]))]
    mu_vFile.write("".join(sList))    
mu_vFile.write("")
mu_vFile.close()

rho_vFile = open("rho_v","w")

for i,T in enumerate(TRange):
    raw_rho_v = str(rho_v[i])[1:-1]
    rho_v1 = raw_rho_v.replace("0","")
    rho_v2 = rho_v1.replace(" ", "")
    rho_v3 = rho_v2.replace(",","")
    raw_pSat = str(pSat[i])[1:-1]
    pSat_l1 = raw_pSat.replace("0","")
    pSat_l2 = pSat_l1.replace(" ", "")
    pSat_l3 = pSat_l2.replace(",","")
    
    sList = ["\t" + rho_v3 + " " + str(Ts[i][j]) + " " +  pSat_l3 + "\n" for j in range(len(Ts[i]))]
    rho_vFile.write("".join(sList))    
rho_vFile.write("")
rho_vFile.close()

rho_lFile = open("rho_l","w")

for i,T in enumerate(TRange):
    raw_rho_l = str(rho_l[i])[1:-1]
    rho_l1 = raw_rho_l.replace("0","")
    rho_l2 = rho_l1.replace(" ", "")
    rho_l3 = rho_l2.replace(",","")
    raw_pSat = str(pSat[i])[1:-1]
    pSat_l1 = raw_pSat.replace("0","")
    pSat_l2 = pSat_l1.replace(" ", "")
    pSat_l3 = pSat_l2.replace(",","")
    
    sList = ["\t" + rho_l3 + " " + str(Ts[i][j]) + " " +  pSat_l3 + "\n" for j in range(len(Ts[i]))]
    rho_lFile.write("".join(sList))    
rho_lFile.write("")
rho_lFile.close()

Cp_lFile = open("Cp_l","w")
Cp_lFile.write("\n")

for i,T in enumerate(TRange):
    raw_Cp_l = str(Cp_l[i])[1:-1]
    Cp_l1 = raw_Cp_l.replace("0","")
    Cp_l2 = Cp_l1.replace(" ", "")
    Cp_l3 = Cp_l2.replace(",","")
    raw_pSat = str(pSat[i])[1:-1]
    pSat_l1 = raw_pSat.replace("0","")
    pSat_l2 = pSat_l1.replace(" ", "")
    pSat_l3 = pSat_l2.replace(",","")
    
    
    Cp_lFile.write("")
    sList = ["\t" + Cp_l3 + " " + str(Ts[i][j]) + " " +  pSat_l3 + "\n" for j in range(len(Ts[i]))]
    Cp_lFile.write("".join(sList))
Cp_lFile.write("")
Cp_lFile.close()

Cp_vFile = open("Cp_v","w")
Cp_vFile.write("\n")

for i,T in enumerate(TRange):
    raw_Cp_v = str(Cp_v[i])[1:-1]
    Cp_v1 = raw_Cp_v.replace("0","")
    Cp_v2 = Cp_v1.replace(" ", "")
    Cp_v3 = Cp_v2.replace(",","")
    raw_pSat = str(pSat[i])[1:-1]
    pSat_l1 = raw_pSat.replace("0","")
    pSat_l2 = pSat_l1.replace(" ", "")
    pSat_l3 = pSat_l2.replace(",","")
    
    Cp_vFile.write("")
    sList = ["\t" + Cp_v3 + " " + str(Ts[i][j]) + " " +  pSat_l3 + "\n" for j in range(len(Ts[i]))]
    Cp_vFile.write("".join(sList))
Cp_vFile.write("")
Cp_vFile.close()

kappa_lFile = open("kappa_l","w")
kappa_lFile.write("\n")

for i,T in enumerate(TRange):
    raw_kappa_l = str(kappa_l[i])[1:-1]
    kappa_l1 = raw_kappa_l.replace("0","")
    kappa_l2 = kappa_l1.replace(" ", "")
    kappa_l3 = kappa_l2.replace(",","")
    raw_pSat = str(pSat[i])[1:-1]
    pSat_l1 = raw_pSat.replace("0","")
    pSat_l2 = pSat_l1.replace(" ", "")
    pSat_l3 = pSat_l2.replace(",","")
    
    kappa_lFile.write("")
    sList = ["\t" + kappa_l3 + " " + str(Ts[i][j]) + " " +  pSat_l3 + "\n" for j in range(len(Ts[i]))]
    kappa_lFile.write("".join(sList))
kappa_lFile.write("")
kappa_lFile.close()

kappa_vFile = open("kappa_v","w")
kappa_vFile.write("\n")

for i,T in enumerate(TRange):
    raw_kappa_v = str(kappa_v[i])[1:-1]
    kappa_v1 = raw_kappa_v.replace("0","")
    kappa_v2 = kappa_v1.replace(" ", "")
    kappa_v3 = kappa_v2.replace(",","")
    raw_pSat = str(pSat[i])[1:-1]
    pSat_l1 = raw_pSat.replace("0","")
    pSat_l2 = pSat_l1.replace(" ", "")
    pSat_l3 = pSat_l2.replace(",","")
    
    kappa_vFile.write("")
    sList = ["\t" + kappa_v3 + " " + str(Ts[i][j]) + " " +  pSat_l3 + "\n" for j in range(len(Ts[i]))]
    kappa_vFile.write("".join(sList))
kappa_vFile.write("")
kappa_vFile.close()

CpMCv_lFile = open("CpMCv_l","w")
CpMCv_lFile.write("\n")

for i,T in enumerate(TRange):
    raw_CpMCv_l = str(CpMCv_l[i])[1:-1]
    CpMCv_l1 = raw_CpMCv_l.replace("0","")
    CpMCv_l2 = CpMCv_l1.replace(" ", "")
    CpMCv_l3 = CpMCv_l2.replace(",","")
    raw_pSat = str(pSat[i])[1:-1]
    pSat_l1 = raw_pSat.replace("0","")
    pSat_l2 = pSat_l1.replace(" ", "")
    pSat_l3 = pSat_l2.replace(",","")
    
    CpMCv_lFile.write("")
    sList = ["\t" + CpMCv_l3 + " " + str(Ts[i][j]) + " " +  pSat_l3 + "\n" for j in range(len(Ts[i]))]
    CpMCv_lFile.write("".join(sList))
CpMCv_lFile.write("")
CpMCv_lFile.close()

CpMCv_vFile = open("CpMCv_v","w")
CpMCv_vFile.write("\n")

for i,T in enumerate(TRange):
    raw_CpMCv_v = str(CpMCv_v[i])[1:-1]
    CpMCv_v1 = raw_CpMCv_v.replace("0","")
    CpMCv_v2 = CpMCv_v1.replace(" ", "")
    CpMCv_v3 = CpMCv_v2.replace(",","")
    raw_pSat = str(pSat[i])[1:-1]
    pSat_l1 = raw_pSat.replace("0","")
    pSat_l2 = pSat_l1.replace(" ", "")
    pSat_l3 = pSat_l2.replace(",","")
    
    CpMCv_vFile.write("")
    sList = ["\t" + CpMCv_v3 + " " + str(Ts[i][j]) + " " +  pSat_l3 + "\n" for j in range(len(Ts[i]))]
    CpMCv_vFile.write("".join(sList))
CpMCv_vFile.write("")
CpMCv_vFile.close()

H_lFile = open("H_l","w")
H_lFile.write("\n")

for i,T in enumerate(TRange):
    raw_H_l = str(H_l[i])[1:-1]
    H_l1 = raw_H_l.replace("0","")
    H_l2 = H_l1.replace(" ", "")
    H_l3 = H_l2.replace(",","")
    raw_pSat = str(pSat[i])[1:-1]
    pSat_l1 = raw_pSat.replace("0","")
    pSat_l2 = pSat_l1.replace(" ", "")
    pSat_l3 = pSat_l2.replace(",","")
     
    H_lFile.write("")
    sList = ["\t" + H_l3 + " " + str(Ts[i][j]) + " " +  pSat_l3 + "\n" for j in range(len(Ts[i]))]
    H_lFile.write("".join(sList))
H_lFile.write("")
H_lFile.close()

H_vFile = open("H_v","w")
H_vFile.write("\n")

for i,T in enumerate(TRange):
    raw_H_v = str(H_v[i])[1:-1]
    H_v1 = raw_H_v.replace("0","")
    H_v2 = H_v1.replace(" ", "")
    H_v3 = H_v2.replace(",","")
    raw_pSat = str(pSat[i])[1:-1]
    pSat_l1 = raw_pSat.replace("0","")
    pSat_l2 = pSat_l1.replace(" ", "")
    pSat_l3 = pSat_l2.replace(",","")
    
    
    H_vFile.write("")
    sList = ["\t" + H_v3 + " " + str(Ts[i][j]) + " " +  pSat_l3 + "\n" for j in range(len(Ts[i]))]
    H_vFile.write("".join(sList))
H_vFile.write("")
H_vFile.close()

E_lFile = open("E_l","w")
E_lFile.write("\n")

for i,T in enumerate(TRange):
    raw_E_l = str(E_l[i])[1:-1]
    E_l1 = raw_E_l.replace("0","")
    E_l2 = E_l1.replace(" ", "")
    E_l3 = E_l2.replace(",","")
    raw_pSat = str(pSat[i])[1:-1]
    pSat_l1 = raw_pSat.replace("0","")
    pSat_l2 = pSat_l1.replace(" ", "")
    pSat_l3 = pSat_l2.replace(",","") 
    
    E_lFile.write("")
    sList = ["\t" + E_l3 + " " + str(Ts[i][j]) + " " +  pSat_l3 + "\n" for j in range(len(Ts[i]))]
    E_lFile.write("".join(sList))
E_lFile.write("")
E_lFile.close()

E_vFile = open("E_v","w")
E_vFile.write("\n")

for i,T in enumerate(TRange):
    raw_E_v = str(E_v[i])[1:-1]
    E_v1 = raw_E_v.replace("0","")
    E_v2 = E_v1.replace(" ", "")
    E_v3 = E_v2.replace(",","")
    raw_pSat = str(pSat[i])[1:-1]
    pSat_l1 = raw_pSat.replace("0","")
    pSat_l2 = pSat_l1.replace(" ", "")
    pSat_l3 = pSat_l2.replace(",","")
    
    E_vFile.write("")
    sList = ["\t" + E_v3 + " " + str(Ts[i][j]) + " " +  pSat_l3 + "\n" for j in range(len(Ts[i]))]
    E_vFile.write("".join(sList))
E_vFile.write("")
E_vFile.close()

S_lFile = open("S_l","w")
S_lFile.write("\n")

for i,T in enumerate(TRange):
    raw_S_l = str(S_l[i])[1:-1]
    S_l1 = raw_S_l.replace("0","")
    S_l2 = S_l1.replace(" ", "")
    S_l3 = S_l2.replace(",","")
    raw_pSat = str(pSat[i])[1:-1]
    pSat_l1 = raw_pSat.replace("0","")
    pSat_l2 = pSat_l1.replace(" ", "")
    pSat_l3 = pSat_l2.replace(",","")
    
    S_lFile.write("")
    sList = ["\t" + S_l3 + " " + str(Ts[i][j]) + " " +  pSat_l3 + "\n" for j in range(len(Ts[i]))]
    S_lFile.write("".join(sList))
S_lFile.write("")
S_lFile.close()

S_vFile = open("S_v","w")
S_vFile.write("\n")

for i,T in enumerate(TRange):
    raw_S_v = str(S_v[i])[1:-1]
    S_v1 = raw_S_v.replace("0","")
    S_v2 = S_v1.replace(" ", "")
    S_v3 = S_v2.replace(",","")
    raw_pSat = str(pSat[i])[1:-1]
    pSat_l1 = raw_pSat.replace("0","")
    pSat_l2 = pSat_l1.replace(" ", "")
    pSat_l3 = pSat_l2.replace(",","")
    
    S_vFile.write("")
    sList = ["\t" + S_v3 + " " + str(Ts[i][j]) + " " +  pSat_l3 + "\n" for j in range(len(Ts[i]))]
    S_vFile.write("".join(sList))
S_vFile.write("")
S_vFile.close()

c_lFile = open("c_l","w")
c_lFile.write("\n")

for i,T in enumerate(TRange):
    raw_c_l = str(c_l[i])[1:-1]
    c_l1 = raw_c_l.replace("0","")
    c_l2 = c_l1.replace(" ", "")
    c_l3 = c_l2.replace(",","")
    raw_pSat = str(pSat[i])[1:-1]
    pSat_l1 = raw_pSat.replace("0","")
    pSat_l2 = pSat_l1.replace(" ", "")
    pSat_l3 = pSat_l2.replace(",","")
    
    c_lFile.write("")
    sList = ["\t" + c_l3 + " " + str(Ts[i][j]) + " " +  pSat_l3 + "\n" for j in range(len(Ts[i]))]
    c_lFile.write("".join(sList))
c_lFile.write("")
c_lFile.close()

c_vFile = open("c_v","w")
c_vFile.write("\n")

for i,T in enumerate(TRange):
    raw_c_v = str(mu_v[i])[1:-1]
    c_v1 = raw_c_v.replace("0","")
    c_v2 = c_v1.replace(" ", "")
    c_v3 = c_v2.replace(",","")
    raw_pSat = str(pSat[i])[1:-1]
    pSat_l1 = raw_pSat.replace("0","")
    pSat_l2 = pSat_l1.replace(" ", "")
    pSat_l3 = pSat_l2.replace(",","")
    
    c_vFile.write("")
    sList = ["\t" + c_v3 + " " + str(Ts[i][j]) + " " +  pSat_l3 + "\n" for j in range(len(Ts[i]))]
    c_vFile.write("".join(sList))
c_vFile.write("")
c_vFile.close()

#Previous dT method to save computational time:
#dT = drho/CP.PropsSI('d(D)/d(P)|T','D',rhoCur,'T',T,fluid_thermo)*CP.PropsSI('d(P)/d(T)|D','D',rhoCur,'T',T,fluid_thermo)

#Previous dP method to save computational time:
#drho/CP.PropsSI('d(D)/d(P)|T','D',rhoPseudoCrit,'T',Tcrit,fluid_thermo)