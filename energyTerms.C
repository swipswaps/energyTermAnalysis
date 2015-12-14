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
    Calculates and writes the enstrophy of the velocity field U.

    The -noWrite option just outputs the max/min values without writing the
    field.

\*---------------------------------------------------------------------------*/

#include "timeSelector.H"
#include "fvc.H"
#include "fvCFD.H"
#include "fixedGradientFvPatchFields.H"
#include <iostream>
#include <fstream>


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{
	Info << " Energy Terms. Version 0.2.1 " << endl; 
	Info << " --------------------------- " << endl;
    Foam::timeSelector::addOptions();
    Foam::argList::addOption
    (
        "meanname",
        "meanname"
    );

    #include "addDictOption.H"

    #include "setRootCase.H"
    #include "createTime.H"
    Foam::instantList timeDirs = Foam::timeSelector::select0(runTime, args);
    #include "createNamedMesh.H"

    #include "readGravitationalAcceleration.H"
	IOdictionary transportProperties
	(
		IOobject
		(
			"transportProperties",
			runTime.constant(),
			mesh,
			IOobject::MUST_READ,
			IOobject::NO_WRITE
		)
	);
    dimensionedScalar beta(transportProperties.lookup("beta"));
    
        

    forAll(timeDirs, timeI)
    {
		
        runTime.setTime(timeDirs[timeI], timeI);
        if (runTime.timeName() == 0) {
			Info << " Skipping " << runTime.timeName() << Foam::endl;
			continue;
		}

        Foam::Info<< "Time = " << runTime.timeName() << Foam::endl;

		
		Info << "\tReading fields " << endl;
		
		volScalarField Tmean(IOobject("Tmean",runTime.timeName(),mesh,IOobject::MUST_READ),mesh); 
		volScalarField T(IOobject("T",runTime.timeName(),mesh,IOobject::MUST_READ),mesh); 
		volScalarField Ttotal(IOobject("Ttotal",runTime.timeName(),mesh,IOobject::NO_READ), T+Tmean);           
		volVectorField Utotal(IOobject("U",runTime.timeName(),mesh,IOobject::MUST_READ), mesh); 

		volVectorField UmeanNew(IOobject("Umovingmean_"+args.option("meanname"),runTime.timeName(),mesh,IOobject::MUST_READ), mesh); 
		volScalarField TmeanNew(IOobject("Tmovingmean_"+args.option("meanname"),runTime.timeName(),mesh,IOobject::MUST_READ), mesh); 
		
		volScalarField p_rgh(IOobject("p_rgh",runTime.timeName(),mesh,IOobject::MUST_READ),mesh); 
		
		volVectorField Uprime = Utotal - UmeanNew; 
		volScalarField Tprime = Ttotal - TmeanNew;
		
		
		// calculating phimean. 
		surfaceScalarField phimean(IOobject("phimean",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE),linearInterpolate(UmeanNew) & mesh.Sf()    );
		surfaceScalarField phi(IOobject("phi",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE),linearInterpolate(Uprime) & mesh.Sf()    );
		
		// get the rho_tag. 
		volScalarField rhok_tag(IOobject("rhok_tag",runTime.timeName(),mesh),- beta*T);
		  
		// get the diffusion coefficient. 
		volTensorField AnisotropicDiffusion(IOobject("AnisotropicDiffusion",runTime.timeName(),	mesh,IOobject::MUST_READ),	mesh    );
		
		/// ------- ///
		Info << "\tCalculating terms" << endl; 
		
		volTensorField UprimeUprime(Uprime*Uprime); 
		
		volTensorField gradUmeanNew(fvc::grad(UmeanNew));
		
		
		volTensorField Mean_To_Kinetic(IOobject("Mean_To_Kinetic_dV",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::AUTO_WRITE),
									   mesh,
									   dimensionedTensor("name", UprimeUprime.dimensions()*	gradUmeanNew.dimensions(), tensor::zero) );
			
		// component wise multlipication. 
		forAll(Mean_To_Kinetic,cellI)
        {
			for (int i =0 ; i < 9 ; i++) {
				Mean_To_Kinetic.internalField()[cellI].component(i) = UprimeUprime.internalField()[cellI].component(i)*gradUmeanNew.internalField()[cellI].component(i)*mesh.V()[cellI]; 
				
			}
		}
		
		
		
		//volScalarField ConversionTerm((IOobject("Kinetic_To_Potential_dV",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::AUTO_WRITE),
			//						   mesh,
				//					   dimensionedScalar("name", Uprime.dimensions()*Tprime.dimensions()*beta.dimensions()*g.dimensions(),pTraits<scalar>::zero) ));
		
		volScalarField ConversionTerm
        (
            IOobject
            (
                "Kinetic_To_Potential_dV",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ
            ),
            Uprime.component(2)*Tprime*g.component(2)*beta
        );
        
        // there is probably a metter way of doing it. 
        forAll(ConversionTerm,cellI) { 
			ConversionTerm[cellI] *= mesh.V()[cellI];
		}
		
		Mean_To_Kinetic.write(); 
        ConversionTerm.write();
        
        //volScalarField ConversionTermdv(ConversionTerm * mesh.V());
        
        //ConversionTerm = (ConversionTerm*mesh.V())();
		
		
		/*
		volVectorField Kinetic_div_U_U(IOobject("Momentum_div_U_U",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::AUTO_WRITE),
										fvc::div(phi, Uprime)        
							   );
							  
		volVectorField Momentum_div_Umean_U(IOobject("Momentum_div_Umean_U",runTime.timeName(),mesh,IOobject::NO_READ),
										fvc::div(phimean, Uprime)        
							   );
							  
		volVectorField Momentum_div_U_Umean(IOobject("Momentum_div_U_Umean",runTime.timeName(),mesh,IOobject::NO_READ),
										fvc::div(phi, Umean)        
							   );
			
		volVectorField Momentum_minus_laplace_D_U(IOobject("Momentum_minus_laplace_D_U",runTime.timeName(),mesh,IOobject::NO_READ),
										-fvc::laplacian(AnisotropicDiffusion,Uprime)
							   );
			
		volVectorField Momentum_minus_grad_p(IOobject("Momentum_minus_grad_p",runTime.timeName(),mesh,IOobject::NO_READ),
										- fvc::grad(p_rgh)
							   );
		
							   
		volVectorField Momentum_g_rhotag(IOobject("Momentum_g_rhotag",runTime.timeName(),mesh,IOobject::NO_READ),
										g*rhok_tag
							   );
							   
		
		volScalarField Energy_div_U_T(IOobject("Energy_div_U_T",runTime.timeName(),mesh,IOobject::NO_READ),
										fvc::div(phi, T)
									 );		
									 
		volScalarField Energy_div_Umean_T(IOobject("Energy_div_Umean_T",runTime.timeName(),mesh,IOobject::NO_READ),
										fvc::div(phimean, T)
									 );		
									 
		volScalarField Energy_div_U_Tmean(IOobject("Energy_div_U_Tmean",runTime.timeName(),mesh,IOobject::NO_READ),
										fvc::div(phi,Tmean)
									 );		
																	 
		volScalarField Energy_minus_laplacian_D_T(IOobject("Energy_minus_laplacian_D_T",runTime.timeName(),mesh,IOobject::NO_READ),
										- fvc::laplacian(AnisotropicDiffusion,T)		
									 );										 								 
	
	
        Momentum_div_U_U.write();
        Momentum_div_Umean_U.write();
        Momentum_div_U_Umean.write();
        Momentum_minus_laplace_D_U.write();
        Momentum_minus_grad_p.write();
        Momentum_g_rhotag.write();
        Energy_div_U_T.write();
        Energy_div_Umean_T.write();
        Energy_div_U_Tmean.write();
        Energy_minus_laplacian_D_T.write();
		
        Foam::Info<< Foam::endl;
		*/													 
								 
	}
	return 0;  
}

// ************************************************************************* //
