/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

Application
    enstrophy

Description
    Calculates and writes the energy terms.
    
    ** Note that we multiply the inner cells with dV and so the values on the boundaries are incorect. 

    The -noWrite option just outputs the max/min values without writing the
    field.

\*---------------------------------------------------------------------------*/
#include "argList.H"
#include "timeSelector.H"
#include "ReadFields.H"
#include "fvCFD.H"

#include "calculateTermBalances.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
/*
 *	Prints the energy calculations in the solver 
 * 
 * */

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    IOdictionary zonesDict
    (
	IOobject
	(
		"zones",
		runTime.constant(),
		mesh,
		IOobject::MUST_READ,
		IOobject::NO_WRITE
    	    )
     );
	
     wordList zones = zonesDict.lookup("zones");
     forAll(zones, zonindx) { 
	calculateTermBalances calculator(mesh,runTime,zones[zonindx]); 
	calculator.printAll();

     }


     return 0;  
}

// ************************************************************************* //
