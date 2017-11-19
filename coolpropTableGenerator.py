"""
Adapted from @author: Luka Denies from TU Delft.
Changelog:
11/12/2017 - Integration of CoolProp
"""

import CoolProp.CoolProp as CP
import numpy as np
import matplotlib.pyplot as plt

#Fluid for rho, cp, cpMcv, h
fluid ='parahydrogen'

#Fluid for thermal conductivity and viscosity models
fluid_sub = 'hydrogen'

T0 = 15 #K
TMax = 100 #K

p0 = 7.5e4 #Pa
pMax = 10e5 #Pa

Tcrit = CP.PropsSI("Tcrit","parahydrogen")

Ts = []
ps = []
pRange = []

rho = []
mu = []
kappa = []
cp = []
h = []
cpMcv = []

i = 0
j = 0

p = p0
T = T0
while p<pMax:
    pRange.append(p)
    TRange = []
    T = T0
    rho.append([0])
    cp.append([0])
    mu.append([0])
    kappa.append([0])
    cpMcv.append([0])
    h.append([0])    
    rho[i][0] = rhoCur = CP.PropsSI('D','T',T,'P',p,fluid)
    cpCur = CP.PropsSI('C','D',rhoCur,'T',T,fluid) 
    cp[i][0] = cpCur
    mu[i][0] = CP.PropsSI('V','D',rhoCur,'T',T,fluid_sub)
    kappa[i][0] = CP.PropsSI('L','D',rhoCur,'T',T,fluid_sub) 
    cpMcv[i][0] = cpCur-CP.PropsSI('O','D',rhoCur,'T',T,fluid) 
    h[i][0] =  CP.PropsSI('H','D',rhoCur,'T',T,fluid)
    TRange.append(T)
    while T<TMax:
        j += 1
        dT = 0.1 #K
        T += dT
        rhoCur = CP.PropsSI('D','T',T,'P',p,fluid)
        rho[i].append(rhoCur)
        cpCur = CP.PropsSI('C','D',rhoCur,'T',T,fluid)
        cp[i].append(cpCur)
        mu[i].append(CP.PropsSI('V','D',rhoCur,'T',T,fluid_sub))
        kappa[i].append(CP.PropsSI('L','D',rhoCur,'T',T,fluid_sub))
        cpMcv[i].append((cpCur-CP.PropsSI('O','D',rhoCur,'T',T,fluid)))
        h[i].append(CP.PropsSI('H','D',rhoCur,'T',T,fluid))
        TRange.append(T)
    i += 1
    ps.append([p]*len(TRange))    
    rhoPseudoCrit = CP.PropsSI('D','T',Tcrit,'P',p,fluid)
    dp = 1000 #Pa 
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

cpFile = open("cp","w")
cpFile.write("( \n")

for i,p in enumerate(pRange):
    cpFile.write("(" + str(p) + "\n(\n")
    sList = ["\t(" + str(Ts[i][j]) + " " + str(cp[i][j]) + ")\n" for j in range(len(Ts[i]))]
    cpFile.write(" ".join(sList))
    cpFile.write(") ) \n")    
cpFile.write(");")
cpFile.close()

kappaFile = open("kappa","w")
kappaFile.write("( \n")

for i,p in enumerate(pRange):
    kappaFile.write("(" + str(p) + "\n(\n")
    sList = ["\t(" + str(Ts[i][j]) + " " + str(kappa[i][j]) + ")\n" for j in range(len(Ts[i]))]
    kappaFile.write(" ".join(sList))
    kappaFile.write(") ) \n")    
kappaFile.write(");")
kappaFile.close()

cpMcvFile = open("cpMcv","w")
cpMcvFile.write("( \n")

for i,p in enumerate(pRange):
    cpMcvFile.write("(" + str(p) + "\n(\n")
    sList = ["\t(" + str(Ts[i][j]) + " " + str(cpMcv[i][j]) + ")\n" for j in range(len(Ts[i]))]
    cpMcvFile.write(" ".join(sList))
    cpMcvFile.write(") ) \n")    
cpMcvFile.write(");")
cpMcvFile.close()

hFile = open("h","w")
hFile.write("( \n")

for i,p in enumerate(pRange):
    hFile.write("(" + str(p) + "\n(\n")
    sList = ["\t(" + str(Ts[i][j]) + " " + str(h[i][j]) + ")\n" for j in range(len(Ts[i]))]
    hFile.write(" ".join(sList))
    hFile.write(") ) \n")    
hFile.write(");")
hFile.close()


#Previous dT method to save computational time:
#dT = drho/CP.PropsSI('d(D)/d(P)|T','D',rhoCur,'T',T,fluid)*CP.PropsSI('d(P)/d(T)|D','D',rhoCur,'T',T,fluid)

#Previous dP method to save computational time:
#drho/CP.PropsSI('d(D)/d(P)|T','D',rhoPseudoCrit,'T',Tcrit,fluid)
