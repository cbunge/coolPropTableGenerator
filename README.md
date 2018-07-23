# tabularCoolProp
CoolProp integrated thermophysical property class for OpenFOAM 5.x

 - To use install CoolProp (https://github.com/CoolProp/CoolProp) 
 - Edit the python file in the coolPropTableGenerator with your fluid and expected T and p ranges.
 - Copy these files into you run/constant folder. Ensure either the enthalpy-based or internal-energy-based tables are implemented in the run/constant file based on your energy type.
 - Copy the /src directory over a fresh copy of the default /src from the original OpenFOAM repository location to your $WM_PROJECT_USER_DIR (this can be found by executing a "echo $WM_PROJECT_USER_DIR" - and the folder that this command identifies may have to be created through "mkdir")
 - Compile the /src/thermophysicalModels/specie folder by executing a "wclean lib" then a "wmake libso"
 - After this is complete compile the /basic folder in the same way.
 - Next compile your solver in the same way as shown in tabularRhoCentralFoam and rhoSimpleFoamTabular. This means altering the title in the /Make/file  and linking the TabularThermophyicalModels and userspecie libs in the EXE_LIBS section. Also note the inclusion of FOAM_USER_LIBBIN.
 - A sample run is shown in the run folder. 

This repository is based off the following work found here. Thank you Luka and Cyril for the help. This is essentially a combination of the extrapolation 2-table design from chriss85 (now continued by Yuusha0) and the interpolation scheme used by ldenies.

For reference, Luka Denies tabulated thermophysicalModels can be found here https://github.com/ldenies/tabulatedProperties and CoolProp here: https://github.com/CoolProp/CoolProp 
