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
    
    ** Note that we multiply the inner cells with dV and so the values on the boundaries are incorect. 

    The -noWrite option just outputs the max/min values without writing the
    field.

\*---------------------------------------------------------------------------*/
#include "argList.H"
#include "timeSelector.H"
#include "ReadFields.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "uniformDimensionedFields.H"
#include "fvCFD.H"

using namespace Foam;

scalar IntegrateVector(const volScalarField & field) { 
		scalar sfield = sum(field).value(); 
		reduce(sfield,sumOp<scalar>());
		return sfield;
}

scalar IntegrateVector(const tmp<volScalarField> & field) { 
		scalar sfield = sum(field).value(); 
		reduce(sfield,sumOp<scalar>());
		return sfield;
}

vector IntegrateVector(const volVectorField & field) { 
		vector sfield = sum(field).value(); 
		reduce(sfield,sumOp<vector>());
		return sfield;
}


vector IntegrateVector(const tmp<volVectorField> & field) { 
		vector sfield = sum(field).value(); 
		reduce(sfield,sumOp<vector>());
		return sfield;
}

vector IntegrateVector(const tmp<volTensorField> & field) { 
		tensor sfield = sum(field).value(); 
		reduce(sfield,sumOp<tensor>());
		return vector(sfield.xx(),sfield.yy(),sfield.zz());
}

tensor IntegrateTensor(const tmp<volTensorField> & field) { 
		tensor sfield = sum(field).value(); 
		reduce(sfield,sumOp<tensor>());
		return sfield;
}

tensor IntegrateTensor(const volTensorField& field) { 
		tensor sfield = sum(field).value(); 
		reduce(sfield,sumOp<tensor>());
		return sfield;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
/*
 *	1.0.0 - 
 * 		Rewriting with OF7  
 * 	0.4.0 - 
 *        Added kinetic energy of the perturbation.  
 * 
 * */

int main(int argc, char *argv[])
{

    Foam::timeSelector::addOptions();
    #include "addRegionOption.H"
    #include "addFunctionObjectOptions.H"

    // Set functionObject post-processing mode
    functionObject::postProcess = true;

    #include "setRootCase.H"

    if (args.optionFound("list"))
    {
        functionObjectList::list();
        return 0;
    }

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
    dimensionedScalar nu(transportProperties.lookup("nu"));
        

    label cellzoneID  =  mesh.cellZones().findZoneID("Work");
    const labelList& cellzonelist =  mesh.cellZones()[cellzoneID];
    scalar tsize = cellzonelist.size() ;
    Info << "zone work has " << tsize << " cells " << endl;
	volScalarField zoneSelector
	(
            IOobject
            (
                    "ZoneSelector",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("zero",dimless,0.0)
	);
    forAll(cellzonelist,cellindx) {
		label currentZoneIndx = cellzonelist[cellindx];
		zoneSelector[currentZoneIndx] = mesh.V()[currentZoneIndx];

    }
    scalar ZoneVolume = sum(zoneSelector).value(); 
    Info << "The zone volume is " << ZoneVolume << endl; 

// =======================================================================
    Foam::surfaceScalarField phi_mean(IOobject("phi_mean",runTime.timeName(),mesh,IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE),
						mesh,dimensionedScalar("tmp",dimLength*dimLength*dimLength/dimTime,0)); 

    Foam::surfaceScalarField phi_tag(IOobject("phi_tag",runTime.timeName(),mesh,IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE),
						mesh,dimensionedScalar("tmp",dimLength*dimLength*dimLength/dimTime,0)); 

    volVectorField U_mean(IOobject("U_mean",runTime.timeName(),mesh,IOobject::NO_READ), mesh,dimensionedVector("a",dimLength/dimTime,vector::zero));
    volScalarField p_mean(IOobject("p_mean",runTime.timeName(),mesh,IOobject::NO_READ), mesh, dimensionedScalar("tmp",dimLength*dimLength/(dimTime*dimTime),0));

    volTensorField mean_to_perb(IOobject("p_mean",runTime.timeName(),mesh,IOobject::NO_READ), mesh, dimensionedTensor("tmp",dimLength*dimLength/(dimTime*dimTime*dimTime),tensor::zero));

    volVectorField U_tag(IOobject("U_tag",runTime.timeName(),mesh,IOobject::NO_READ), mesh,dimensionedVector("a",dimLength/dimTime,vector::zero));
    volScalarField p_tag(IOobject("p_tag",runTime.timeName(),mesh,IOobject::NO_READ), mesh, dimensionedScalar("tmp",dimLength*dimLength/(dimTime*dimTime),0));

    volTensorField reynoldsU(IOobject("reynoldsU",runTime.timeName(),mesh,IOobject::NO_READ), mesh,dimensionedTensor("tmp",dimLength*dimLength/(dimTime*dimTime),tensor::zero));

    volVectorField U_mean_dt(IOobject("U_mean_dt",runTime.timeName(),mesh,IOobject::NO_READ), mesh,dimensionedVector("a",dimLength/dimTime,vector::zero));

// =======================================================================

    scalar DT = timeDirs[timeDirs.size()-1].value() - timeDirs[0].value();
    dimensionedScalar DeltaT("dt",dimTime,DT);

// ============================ Calculate the mean =======================
    Info << " Calculating Mean" << endl;

    dimensionedScalar mean_Energy_All("mm", dimLength*dimLength/(dimTime*dimTime*dimTime),0);
    dimensionedScalar mean_EnergyDT_All("mm", dimLength*dimLength/(dimTime*dimTime*dimTime),0);
    dimensionedScalar mean_EnergyAdvection_All("mm", dimLength*dimLength/(dimTime*dimTime*dimTime),0);
    dimensionedScalar mean_EnergyDiffusion_All("mm", dimLength*dimLength/(dimTime*dimTime*dimTime),0);
    dimensionedScalar mean_EnergyPressureConvection_All("mm", dimLength*dimLength/(dimTime*dimTime*dimTime),0);
    dimensionedScalar mean_EnergyT_All("mm", dimLength*dimLength/(dimTime*dimTime*dimTime),0);

    dimensionedScalar mean_Energy_Zone("mm", dimLength*dimLength/(dimTime*dimTime*dimTime),0);
    dimensionedScalar mean_EnergyDT_Zone("mm", dimLength*dimLength/(dimTime*dimTime*dimTime),0);
    dimensionedScalar mean_EnergyAdvection_Zone("mm", dimLength*dimLength/(dimTime*dimTime*dimTime),0);
    dimensionedScalar mean_EnergyDiffusion_Zone("mm", dimLength*dimLength/(dimTime*dimTime*dimTime),0);
    dimensionedScalar mean_EnergyPressureConvection_Zone("mm", dimLength*dimLength/(dimTime*dimTime*dimTime),0);


    forAll(timeDirs, timeI)
    {
	Info << timeDirs[timeI] << endl;
        runTime.setTime(timeDirs[timeI], timeI);
	surfaceScalarField phi(IOobject("phi",runTime.timeName(),mesh,IOobject::MUST_READ,IOobject::NO_WRITE),mesh);
	volVectorField U(IOobject("U",runTime.timeName(),mesh,IOobject::MUST_READ), mesh);
 	volScalarField p(IOobject("p",runTime.timeName(),mesh,IOobject::MUST_READ), mesh);

	if (timeI == 0) { 
		U_mean_dt = -U;
		continue;
	}

	word prevTime= timeDirs[timeI-1].name();
	dimensionedScalar dt("dt",dimTime,timeDirs[timeI].value() - timeDirs[timeI-1].value()); 
	volVectorField Uold(IOobject("U",prevTime,mesh,IOobject::MUST_READ), mesh);
	//Info << "dt " << dt << endl;

	phi_mean += phi;
	U_mean   += U; 
	p_mean   += p;

	
	if (timeI == timeDirs.size()-1) {
		U_mean_dt += U;
	}


	mean_Energy_All 			+= (U&((U-Uold)/dt+fvc::div(phi,U)-fvc::laplacian(nu,U)+fvc::grad(p)))->weightedAverage(mesh.V());
	mean_EnergyDT_All			+= (U&(U-Uold)/dt)->weightedAverage(mesh.V());
	mean_EnergyAdvection_All		+= (U&fvc::div(phi,U))->weightedAverage(mesh.V()); 
	mean_EnergyDiffusion_All		+= (U&fvc::laplacian(nu,U))->weightedAverage(mesh.V());
	mean_EnergyPressureConvection_All	+= (U&fvc::grad(p))->weightedAverage(mesh.V());


/*
	Info << "-------Zone---------" << endl;
	Info << "1" << endl;
	Info << dt << endl;
	Info << (U&(U-Uold)/dt)->dimensions() << endl;
	mean_EnergyDT_Zone += IntegrateVector( (U&(U-Uold)/dt)*zoneSelector);
	
	Info << "2" << endl;
	mean_EnergyAdvection_Zone += IntegrateVector( (U&fvc::div(phi,U))*zoneSelector);
	Info << "3" << endl;	
	mean_EnergyDiffusion_Zone += IntegrateVector( (U&fvc::laplacian(nu,U))*zoneSelector);
	Info << "4" << endl;
	mean_EnergyPressureConvection_Zone += IntegrateVector( (U&fvc::grad(p))*zoneSelector); 
	Info << "5" << endl;
	mean_Energy_Zone += IntegrateVector( (U&((U-Uold)/dt+fvc::div(phi,U)-fvc::laplacian(nu,U)+fvc::grad(p)))*zoneSelector);

*/

    }
    phi_mean  /= timeDirs.size();
    U_mean    /= timeDirs.size();
    p_mean    /= timeDirs.size();
    U_mean_dt /= DeltaT;


    mean_Energy_All/= timeDirs.size();
    mean_EnergyDT_All/= timeDirs.size();
    mean_EnergyAdvection_All/= timeDirs.size();
    mean_EnergyDiffusion_All/= timeDirs.size();
    mean_EnergyPressureConvection_All/= timeDirs.size();

    mean_Energy_Zone		/= timeDirs.size();
    mean_EnergyDT_Zone		/= timeDirs.size();
    mean_EnergyAdvection_Zone	/= timeDirs.size();
    mean_EnergyDiffusion_Zone	/= timeDirs.size();
    mean_EnergyPressureConvection_Zone/= timeDirs.size();



	//Info << " First calc " << mean_to_perb << " || " << (reynoldsU&&fvc::grad(U_mean))->weightedAverage(mesh.V()) << endl;

// =======================================================================

	Info << " =========================== Mean balances ================= " << endl; 
	Info << "Momentum balance" << endl; 

	Info << "-------All domain---------" << endl;
	Info << "DT " 	     	<< U_mean_dt.weightedAverage(mesh.V()) << endl;
	Info << "Advection " 	<< (fvc::div(phi_mean,U_mean)+fvc::div(reynoldsU) )->weightedAverage(mesh.V()) << endl;
	Info << "Diffusion " 	<< fvc::laplacian(nu,U_mean)->weightedAverage(mesh.V()) << endl;
	Info << "Pressure  " 	<< fvc::grad(p_mean)->weightedAverage(mesh.V()) << endl;
	Info << "Total " 	<< (U_mean_dt + fvc::div(phi_mean,U_mean)+fvc::div(reynoldsU) - fvc::laplacian(nu,U_mean) + fvc::grad(p_mean) )->weightedAverage(mesh.V()) << endl;

	Info << "-------Zone---------" << endl;
	Info << "DT " 	     << IntegrateVector( U_mean_dt*zoneSelector) << endl;
	Info << "Advection " << IntegrateVector( (fvc::div(phi_mean,U_mean)+fvc::div(reynoldsU))*zoneSelector) << endl;
	Info << "Diffusion " << IntegrateVector( fvc::laplacian(nu,U_mean)*zoneSelector) << endl;
	Info << "Pressure  " << IntegrateVector( fvc::grad(p_mean)*zoneSelector) << endl;
	Info << "Total " << IntegrateVector( (U_mean_dt + fvc::div(phi_mean,U_mean)+fvc::div(reynoldsU) - fvc::laplacian(nu,U_mean) + fvc::grad(p_mean) )*zoneSelector) << endl;

	Info << "Energy balance" << endl; 
	Info << "-------All domain---------" << endl;
	Info << "DT " 	     << (U_mean&U_mean_dt)->weightedAverage(mesh.V()) << endl;
	Info << "Advection " << (U_mean&(fvc::div(phi_mean,U_mean)+fvc::div(reynoldsU) ))->weightedAverage(mesh.V()) << endl;
	Info << "Diffusion " << (U_mean&fvc::laplacian(nu,U_mean))->weightedAverage(mesh.V()) << endl;
	Info << "Pressure  " << (U_mean&fvc::grad(p_mean))->weightedAverage(mesh.V()) << endl;
	Info << "Total " << (U_mean&(U_mean_dt + fvc::div(phi_mean,U_mean)+fvc::div(reynoldsU) - fvc::laplacian(nu,U_mean) + fvc::grad(p_mean) ))->weightedAverage(mesh.V()) << endl;


	Info << endl << endl;

// =======================================================================
    dimensionedScalar MeanPerturbation_Energy("mm", dimLength*dimLength/(dimTime*dimTime*dimTime),0);
    dimensionedScalar meanPerturbation_EnergyDT("mm", dimLength*dimLength/(dimTime*dimTime*dimTime),0);
    dimensionedScalar meanPerturbation_EnergyAdvection("mm", dimLength*dimLength/(dimTime*dimTime*dimTime),0);
    dimensionedScalar meanPerturbation_EnergyDiffusion("mm", dimLength*dimLength/(dimTime*dimTime*dimTime),0);
    dimensionedScalar meanPerturbation_EnergyPressureConvection("mm", dimLength*dimLength/(dimTime*dimTime*dimTime),0);

    forAll(timeDirs, timeI)
    {
	if (timeI ==0){
		continue;
	}
		
        runTime.setTime(timeDirs[timeI], timeI);
        if (runTime.time() == 0) {
			Info << " Skipping " << runTime.timeName() << Foam::endl;
			continue;
		}

        Foam::Info<< "Time = " << runTime.timeName() << Foam::endl;

	surfaceScalarField phi(IOobject("phi",runTime.timeName(),mesh,IOobject::MUST_READ,IOobject::NO_WRITE),mesh);
	volVectorField U(IOobject("U",runTime.timeName(),mesh,IOobject::MUST_READ), mesh);
 	volScalarField p(IOobject("p",runTime.timeName(),mesh,IOobject::MUST_READ), mesh);

	word prevTime= timeDirs[timeI-1].name();
	dimensionedScalar dt("dt",dimTime,timeDirs[timeI].value() - timeDirs[timeI-1].value()); 
	volVectorField Uold(IOobject("U",prevTime,mesh,IOobject::MUST_READ), mesh);
	Info << "dt " << dt << endl;
	Info << endl;
//	----------------------------------------------------- Reynolds
	surfaceScalarField phi_tag = phi-phi_mean;
	volVectorField tagU 	= U - U_mean;
	volVectorField tagU_old = Uold - U_mean;
	volScalarField tagp 	= p - p_mean;
	
	volTensorField reynoldsU = tagU*tagU;

//	------------------------------------------------------------ Utag terms. 
	volVectorField U_tag_convection  = fvc::div(phi,tagU);
	volVectorField U_convection_Utag = fvc::div(phi_tag,U_mean);
	volVectorField U_tag_laplacian   = fvc::laplacian(nu,tagU);
	volVectorField p_tag_grad        = fvc::grad(tagp); 
	volVectorField U_tag_dt_diff 	 = (U - Uold)/dt- U_mean_dt.weightedAverage(mesh.V());
	
	reynoldsU += tagU*tagU;
/*	volTensorField reynoldsUTS(IOobject("reynoldsU",runTime.timeName(),mesh,IOobject::NO_READ), tagU*tagU);
	reynoldsUTS.write();
*/
	mean_to_perb += (reynoldsU & fvc::grad(U_mean) ); 
//	--------------------------------------------------------------------	
	meanPerturbation_EnergyDT		  += (U&U_tag_dt_diff)->weightedAverage(mesh.V());;
	meanPerturbation_EnergyPressureConvection += (fvc::div(phi,0.5*(U&U)))->weightedAverage(mesh.V());
	meanPerturbation_EnergyAdvection	  += (U&fvc::div(phi,U))->weightedAverage(mesh.V());
//	----------------------------------------------------- Prints
	Info << "Continuity" << endl; 
	Info << "----------" << endl;
	Info << fvc::div(phi)->weightedAverage(mesh.V()) << endl;

	Info << endl;
	Info << "Momentum balance" << endl; 
	Info << "-------All domain---------" << endl;

	Info << "DT " 	     << ((U-Uold)/dt)->weightedAverage(mesh.V()) << endl;
	Info << "DT diff "   << ((U-Uold)/dt)->weightedAverage(mesh.V()) - U_mean_dt.weightedAverage(mesh.V()) << endl;
	Info << "Advection " << fvc::div(phi,U)->weightedAverage(mesh.V()) << endl;
	Info << "Diffusion " << fvc::laplacian(nu,U)->weightedAverage(mesh.V()) << endl;
	Info << "Pressure  " << fvc::grad(p)->weightedAverage(mesh.V()) << endl;
	Info << "Total " << ((U-Uold)/dt+fvc::div(phi,U)-fvc::laplacian(nu,U)+fvc::grad(p))->weightedAverage(mesh.V()) << endl;


	Info << endl;
	Info << "-------Zone---------" << endl;
	Info << "DT " 	     << IntegrateVector( ((U-Uold)/dt)*zoneSelector) << endl;
	Info << "Advection " << IntegrateVector( fvc::div(phi,U)*zoneSelector) << endl;
	Info << "Diffusion " << IntegrateVector( fvc::laplacian(nu,U)*zoneSelector)<< endl;
	Info << "Pressure  " << IntegrateVector( fvc::grad(p)*zoneSelector) << endl;
	Info << "Total "     << IntegrateVector( ((U-Uold)/dt+fvc::div(phi,U)-fvc::laplacian(nu,U)+fvc::grad(p))*zoneSelector) << endl;

	Info << endl << " -=-=-=-=-=-=-=-=-=-= Perturbation "  << endl;
	Info << "-------All domain---------" << endl;
	Info << "DT " << U_tag_dt_diff.weightedAverage(mesh.V()) << endl;
	Info << "div(U,tagU) " << U_tag_convection.weightedAverage(mesh.V()) << endl;
	Info << "div(Umean,tagU) " << U_convection_Utag.weightedAverage(mesh.V()) << endl;
	Info << "-overline(div(tagU,tagU)) " << fvc::div(reynoldsU)->weightedAverage(mesh.V()) << endl;
	Info << "Advection " << (U_tag_convection + U_convection_Utag - fvc::div(reynoldsU) )->weightedAverage(mesh.V()) << endl;
	Info << "Diffusion  " << U_tag_laplacian.weightedAverage(mesh.V()) << endl;
	Info << "Pressure " << p_tag_grad.weightedAverage(mesh.V()) << endl;
	Info << "Total " << (U_tag_dt_diff+U_tag_convection+U_convection_Utag-fvc::div(reynoldsU)-U_tag_laplacian+p_tag_grad)->weightedAverage(mesh.V()) << endl;
	Info << endl;
	Info << "-------Zone---------" << endl;
	Info << "DT " 				<< IntegrateVector(U_tag_dt_diff*zoneSelector) << endl;
	Info << "div(U,tagU) " 			<< IntegrateVector(U_tag_convection*zoneSelector) << endl;
	Info << "div(Umean,tagU) " 		<< IntegrateVector(U_convection_Utag*zoneSelector) << endl;
	Info << "-overline(div(tagU,tagU)) " 	<< IntegrateVector(fvc::div(reynoldsU)*zoneSelector) << endl;
	Info << "Advection " 			<< IntegrateVector((U_tag_convection + U_convection_Utag - fvc::div(reynoldsU) )*zoneSelector) << endl;
	Info << "Diffusion  " 			<< IntegrateVector(U_tag_laplacian*zoneSelector) << endl;
	Info << "Pressure " 			<< IntegrateVector(p_tag_grad*zoneSelector) << endl;
	Info << "Total " 			<< IntegrateVector((U_tag_dt_diff+U_tag_convection+U_convection_Utag-fvc::div(reynoldsU)-U_tag_laplacian+p_tag_grad)*zoneSelector) << endl;

	Info <<endl;
	Info << "Energy balance" << endl; 
	Info << "-------All domain---------" << endl;
	Info << "DT " 	     << (U&(U-Uold)/dt)->weightedAverage(mesh.V()) << endl;
	Info << "Advection " << (U&fvc::div(phi,U))->weightedAverage(mesh.V()) << endl;
	Info << "Diffusion " << (U&fvc::laplacian(nu,U))->weightedAverage(mesh.V()) << endl;
	Info << "Pressure  " << (U&fvc::grad(p))->weightedAverage(mesh.V()) << endl;
	Info << "Pressure  2" << (fvc::div(phi,p))->weightedAverage(mesh.V()) << endl;
	Info << "Total " << (U&((U-Uold)/dt+fvc::div(phi,U)-fvc::laplacian(nu,U)+fvc::grad(p)))->weightedAverage(mesh.V()) << endl;


	Info << "-------Zone---------" << endl;
	Info << "DT " 	     << IntegrateVector( (U&(U-Uold)/dt)*zoneSelector) << endl;
	Info << "Advection " << IntegrateVector( (U&fvc::div(phi,U))*zoneSelector) << endl;
	Info << "Diffusion " << IntegrateVector( (U&fvc::laplacian(nu,U))*zoneSelector)<< endl;
	Info << "Pressure  " << IntegrateVector( (U&fvc::grad(p))*zoneSelector) << endl;
	Info << "Total "     << IntegrateVector( (U&((U-Uold)/dt+fvc::div(phi,U)-fvc::laplacian(nu,U)+fvc::grad(p)))*zoneSelector) << endl;


	Info << endl << " -=-=-=-=-=-=-=-=-=-= Perturbation "  << endl;
	Info << "-------All domain---------" << endl;
	Info << "DT " << (U&U_tag_dt_diff)->weightedAverage(mesh.V()) << endl;
	Info << "div(U,tagU) " << (U&U_tag_convection)->weightedAverage(mesh.V()) << endl;
	Info << "div(Umean,tagU) " << (U&U_convection_Utag)->weightedAverage(mesh.V()) << endl;
	Info << "-overline(div(tagU,tagU)) " << (U&fvc::div(reynoldsU))->weightedAverage(mesh.V()) << endl;
	Info << "Advection " << (U&(U_tag_convection + U_convection_Utag - fvc::div(reynoldsU) ))->weightedAverage(mesh.V()) << endl;
	Info << "Diffusion  " << (U&U_tag_laplacian)->weightedAverage(mesh.V()) << endl;
	Info << "Pressure " << (U&p_tag_grad)->weightedAverage(mesh.V()) << endl;
	Info << "Total " << (U&(U_tag_dt_diff+U_tag_convection+U_convection_Utag-fvc::div(reynoldsU)-U_tag_laplacian+p_tag_grad))->weightedAverage(mesh.V()) << endl;

	Info <<endl << endl;

	}

    	reynoldsU /= timeDirs.size();
	mean_to_perb /= timeDirs.size();

	Info << " Final perturbation statistics " << endl;


	Info << "Mean energy u'du'/dt " << 	mean_EnergyDT_All/timeDirs.size()-(U_mean&U_mean_dt)->weightedAverage(mesh.V())  << endl;
	Info << "Mean energy advection+[transger to perb] " <<	mean_EnergyAdvection_All/timeDirs.size()-(U_mean&(fvc::div(phi_mean,U_mean) ))->weightedAverage(mesh.V()) << endl; //+
	Info << "Mean energy tranfer to perb" << 	mean_to_perb.weightedAverage(mesh.V()) << endl;
	Info << "Mean energy diffusion " <<	mean_EnergyDiffusion_All/timeDirs.size()-(U_mean&fvc::laplacian(nu,U_mean))->weightedAverage(mesh.V()) << endl;
	Info << "Mean energy pressure " <<	mean_EnergyPressureConvection_All/timeDirs.size()- (U_mean&fvc::grad(p_mean))->weightedAverage(mesh.V())<< endl;

	
	return 0;  
}

// ************************************************************************* //


/*		
		Info << "\tReading fields " << endl;
		
		volScalarField Tmean(IOobject("Tmean",runTime.timeName(),mesh,IOobject::MUST_READ),mesh); 
		volScalarField T(IOobject("T",runTime.timeName(),mesh,IOobject::MUST_READ),mesh); 
		volScalarField Ttotal(IOobject("Ttotal",runTime.timeName(),mesh,IOobject::NO_READ), T+Tmean);           
		volVectorField Utotal(IOobject("U",runTime.timeName(),mesh,IOobject::MUST_READ), mesh); 

		volVectorField UmeanNew(IOobject("Umovingmean_"+args.option("meanname"),runTime.timeName(),mesh,IOobject::MUST_READ), mesh); 
		volScalarField TmeanNew(IOobject("Tmovingmean_"+args.option("meanname"),runTime.timeName(),mesh,IOobject::MUST_READ), mesh); 
		
		volScalarField p_rgh(IOobject("p_rgh",runTime.timeName(),mesh,IOobject::MUST_READ),mesh); 
		
		volVectorField Uprime(IOobject("Uprime",runTime.timeName(),	mesh,IOobject::NO_READ,IOobject::AUTO_WRITE),Utotal - UmeanNew); 
		volScalarField Tprime = T      - TmeanNew;
		
		
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
		
		
		volTensorField UgradVec1(IOobject("UgradD",runTime.timeName(),	mesh,IOobject::NO_READ,IOobject::AUTO_WRITE),AnisotropicDiffusion&fvc::grad(Uprime));
		volTensorField UgradVec2(IOobject("Ugrad",runTime.timeName(),	mesh,IOobject::NO_READ,IOobject::AUTO_WRITE),	fvc::grad(Uprime));
		
		
		volScalarField KineticDiffusion(IOobject("KineticEnergyDiffusion_dV",
												 runTime.timeName(),mesh,IOobject::NO_READ,IOobject::AUTO_WRITE),
												 UgradVec1&&UgradVec2);
		
		//volScalarField V(IOobject("Vmesh",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::AUTO_WRITE),mesh,dimensionedScalar("name", mesh.V().dimensions(), scalar(0)) );
		
		volTensorField Mean_To_Kinetic(IOobject("Mean_To_Kinetic_dV",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::AUTO_WRITE),
									   mesh,
									   dimensionedTensor("name", UprimeUprime.dimensions()*	gradUmeanNew.dimensions(), tensor::zero) );
		
		volScalarField ConversionTerm
        (
            IOobject
            (
                "Kinetic_To_Potential_dV",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ
            ),
            -Uprime.component(2)*Tprime*g.component(2)*beta // remember that the g is with + sign. 
        );

        volScalarField PrimeKineticEnergy
        (
            IOobject
            (
                "PrimeKineticEnergy_dV",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ
            ),
            0.5*magSqr(Uprime)
        );
        
        
        
        // there is probably a metter way of doing it. 
        // !!! Note that we multiply the inner cells with dV and so the values on the boundaries are incorect. !!!
        forAll(ConversionTerm,cellI) { 
			
			for (int i =0 ; i < 9 ; i++) {
				Mean_To_Kinetic.internalField()[cellI].component(i) = UprimeUprime.internalField()[cellI].component(i)*gradUmeanNew.internalField()[cellI].component(i)*mesh.V()[cellI]; 
			}
			
			// there is probably a metter way of doing it. 
			ConversionTerm[cellI] *= mesh.V()[cellI];
			
			// there is probably a metter way of doing it. 
			KineticDiffusion[cellI] *= mesh.V()[cellI];
			
			PrimeKineticEnergy[cellI] *= mesh.V()[cellI];
			
			//V[cellI] = mesh.V()[cellI];
		}
		
		
		PrimeKineticEnergy.write();
		Mean_To_Kinetic.write(); 
        ConversionTerm.write();
        KineticDiffusion.write();
        
		
		
		
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


DT (U_mean&U_mean_dt).weightedAverage(weights) [0 2 -3 0 0 0 0] 5.96823
Advection (U_mean&(convection(phi_mean,U_mean)+div(reynoldsU))).weightedAverage(weights) [0 2 -3 0 0 0 0] -7.00016
Diffusion (U_mean&laplacian(nu,U_mean)).weightedAverage(weights) [0 2 -3 0 0 0 0] -0.0102145
Pressure  (U_mean&grad(p_mean)).weightedAverage(weights) [0 2 -3 0 0 0 0] 0.107472
Total (U_mean&((((U_mean_dt+convection(phi_mean,U_mean))+div(reynoldsU))-laplacian(nu,U_mean))+grad(p_mean))).weightedAverage(weights) [0 2 -3 0 0 0 0] -0.914243



Energy balance 2
DT mm [0 2 -3 0 0 0 0] 6.42525
Advection mm [0 2 -3 0 0 0 0] -6.86761 || mm [0 2 -3 0 0 0 0] -6.54538
Diffusion mm [0 2 -3 0 0 0 0] -0.0103222
Pressure mm [0 2 -3 0 0 0 0] 0.122995 || mm [0 2 -3 0 0 0 0] 1.6023e-15


*/
								 
