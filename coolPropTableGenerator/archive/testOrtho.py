"""
Carl Bunge
Washington State University
June 2018

Adapted from @author: Luka Denies from TU Delft.

Changelog:
11/2017 - Integration of CoolProp
06/2018 - Update to OpenFOAM-5.x (Mass-based thermodynamics (for example: cpMcv to CpMCv))
03/2019 - Update to include orthohydrogen extrapolated thermal conductivity
"""

import CoolProp.CoolProp as CP
import numpy as np
import matplotlib.pyplot as plt

#Fluid for thermodynamic properties (rho, Cp, CpMcv, H, E, S, c, thermal conductivity)
fluid_thermo ='orthohydrogen'

print CP.PropsSI('T','P',101325,'Q',0,'OrthoHydrogen')

print CP.PropsSI('Dmolar','P',101325,'Q',0,'OrthoHydrogen')

print CP.PropsSI('Hmolar','P',101325,'Q',0,'OrthoHydrogen')

print CP.PropsSI('Smolar','P',101325,'Q',0,'OrthoHydrogen')


print CP.PropsSI('T','P',101325,'Q',0,'parahydrogen')

print CP.PropsSI('Dmolar','P',101325,'Q',0,'parahydrogen')

print CP.PropsSI('Hmolar','P',101325,'Q',0,'parahydrogen')

print CP.PropsSI('Smolar','P',101325,'Q',0,'parahydrogen')



#Previous dT method to save computational time:
#dT = drho/CP.PropsSI('d(D)/d(P)|T','D',rhoCur,'T',T,fluid_thermo)*CP.PropsSI('d(P)/d(T)|D','D',rhoCur,'T',T,fluid_thermo)

#Previous dP method to save computational time:
#drho/CP.PropsSI('d(D)/d(P)|T','D',rhoPseudoCrit,'T',Tcrit,fluid_thermo)
