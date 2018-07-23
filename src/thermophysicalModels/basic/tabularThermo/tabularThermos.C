/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 Yuusha, tilasoldo and cbunge
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of tilasoldo and Yuusha contribution to OpenFOAM.
    It is based on chriss85 contribution for OpenFOAM 2.3.x.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "psiThermo.H"
#include "rhoThermo.H"
#include "makeThermo.H"

#include "specie.H"
#include "perfectGas.H"
#include "tabularEOS.H"
#include "hConstThermo.H"
#include "hPolynomialThermo.H"
#include "hTabularThermo.H"
#include "eConstThermo.H"
#include "janafThermo.H"
#include "sensibleEnthalpy.H"
#include "sensibleInternalEnergy.H"
#include "thermo.H"

#include "constTransport.H"
#include "sutherlandTransport.H"
#include "polynomialTransport.H"
#include "tabularTransport.H"

#include "hePsiThermo.H"
#include "heRhoThermo.H"
#include "heTabularThermo.H"
#include "pureMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/* * * * * * * * * * * * * * Internal-energy-based * * * * * * * * * * * * * */
// This uses two tables for e(p, T) and T (e, p) NOTE: Ensure the correct table (internalEnergy or Enthalpy) is being referened in run/constant/thermophysicalProperties file.
makeThermo
(
    psiThermo,
    heTabularThermo,
    pureMixture,
    tabularTransport,
    sensibleInternalEnergy,
    hTabularThermo,
    tabularEOS,
    specie
);
 
/* * * * * * * * * * * * * * Enthalpy-based * * * * * * * * * * * * */
// This uses two tables for h(p, T) and T(h, p) NOTE: Ensure the correct table (internalEnergy or Enthalpy) is being referened in run/constant/thermophysicalProperties file.
makeThermo
(
    psiThermo,
    heTabularThermo,
    pureMixture,
    tabularTransport,
    sensibleEnthalpy,
    hTabularThermo,
    tabularEOS,
    specie
);

   
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
