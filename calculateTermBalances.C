#include "calculateTermBalances.H"

#define getFields(name,type)\
	vol##type##Field total(\
            IOobject		\
            (			\
                    "total_energy_#name",		\
                    word(runTime.endTime().value()),	\
                    mesh,			\
                    IOobject::MUST_READ,		\
                    IOobject::NO_WRITE	\
            ),\
	    mesh\
	); \
	vol##type##Field mean(\
            IOobject		\
            (			\
                    "mean_energy_#name",		\
                    word(runTime.startTime().value()),	\
                    mesh,			\
                    IOobject::MUST_READ,		\
                    IOobject::NO_WRITE	\
            ),\
	    mesh\
	); \
	vol##type##Field perb(\
            IOobject		\
            (			\
		    "perb_energy_#name",		\
                    word(runTime.endTime().value()),	\
                    mesh,			\
                    IOobject::MUST_READ,		\
                    IOobject::NO_WRITE	\
            ),\
	    mesh\
	); 


#define getTensor(name) \
	volTensorField name(\
            IOobject		\
            (			\
                    #name,		\
                    word(runTime.endTime().value()),	\
                    mesh,			\
                    IOobject::MUST_READ,		\
                    IOobject::NO_WRITE	\
            ),\
	    mesh\
	); \


calculateTermBalances::calculateTermBalances(fvMesh& imesh,Time& irunTime, word zoneName) : 
		mesh(imesh),
		runTime(irunTime),
		zoneSelector(
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
		)
{ 

    if (zoneName == "") { 
	    forAll(mesh.V(),cellindx) {
			zoneSelector[cellindx] = mesh.V()[cellindx];
	    }

	    ZoneVolume = sum(zoneSelector).value(); 
	    reduce(ZoneVolume,sumOp<scalar>());
    } else { 

	    label cellzoneID  =  mesh.cellZones().findZoneID(zoneName);
	    const labelList& cellzonelist =  mesh.cellZones()[cellzoneID];
	    //scalar tsize = cellzonelist.size();
		 
	    forAll(cellzonelist,cellindx) {
			label currentZoneIndx = cellzonelist[cellindx];
			zoneSelector[currentZoneIndx] = mesh.V()[currentZoneIndx];
	    }

	    ZoneVolume = sum(zoneSelector).value(); 
	    reduce(ZoneVolume,sumOp<scalar>());
    }

}


void calculateTermBalances::printTemporal() { 
	scalar timeSpan = (runTime.endTime() - runTime.startTime()).value();
	getFields("dUdt",Scalar)

	scalar IntegratedTotal = Integrate(total)/timeSpan;
	scalar IntegratedMean = Integrate(mean)/timeSpan;
	scalar IntegratedPerb = Integrate(perb)/timeSpan;

	Info << "Temporal terms" << endl;
	Info << "--------------" << endl;
	Info << "\t total = mean + perb => " << IntegratedTotal << " =  " << IntegratedMean << " + " << IntegratedPerb  << endl;
	Info << "-----------------------------------" <<endl << endl;	
} //.. printTemporal. 

void calculateTermBalances::printDiffusion() { 
	scalar timeSpan = (runTime.endTime() - runTime.startTime()).value();
	getFields("diffusion",Tensor)

	tensor IntegratedTotal = Integrate(total)/timeSpan;
	tensor IntegratedMean = Integrate(mean)/timeSpan;
	tensor IntegratedPerb = Integrate(perb)/timeSpan;


	
	Info << "Total Diffusion term (before partial integration)" << endl;
	Info << "-------------------------------------------------" << endl;
	Info << "total = mean + perb" << endl;

	word component;
	word dir1;
	word dir2;
	
	label c=0;
	for (int i=0 ; i < 3 ; i++ ) {
		switch (i) {
			case 0: 
				dir1 = "x";
				break;
			case 1: 
				dir1 = "y";
				break;
			case 2: 
				dir1 = "z";
				break;
		};
		for (int j=0 ; j < 3 ; j++,c++ ) {
			switch (j) {
				case 0: 
					component="u";
					dir2 = "x";
					break;
				case 1: 
					component="v";
					dir2 = "y";
					break;
				case 2: 
					component="w";
					dir2 = "z";
					break;
			};

			Info << "\t\t d" << component << "/d"<< dir1 << dir1 << 
					" total " << IntegratedTotal.component(c) << " = " 
						  << IntegratedMean.component(c)  << " + " << IntegratedPerb.component(c)  << endl;
					
		} //..for 
	} //.. for 
	Info << "-----------------------------------" <<endl << endl;	
} //.. printDiffusion 

void calculateTermBalances::printPotential() { 

	scalar timeSpan = (runTime.endTime() - runTime.startTime()).value();
	getFields("potential",Scalar)

	scalar IntegratedTotal = Integrate(total)/timeSpan;
	scalar IntegratedMean = Integrate(mean)/timeSpan;
	scalar IntegratedPerb = Integrate(perb)/timeSpan;

	Info << "Potential" << endl;
	Info << "---------" << endl; 
	Info << "total = mean + perb => " << IntegratedTotal  << " = " <<  IntegratedMean << "+ " <<  IntegratedPerb << endl;
	Info << "-----------------------------------" <<endl << endl;	

}

void calculateTermBalances::printKUsqr() { 

	scalar timeSpan = (runTime.endTime() - runTime.startTime()).value();
	volTensorField mean(
            IOobject		
            (			
                    "mean_KgradUsqr",		
                    word(runTime.startTime().value()),	
                    mesh,			
                    IOobject::MUST_READ,		
                    IOobject::NO_WRITE	
            ),
	    mesh
	); 

	volTensorField perb(
            IOobject		
            (			
                    "perb_KgradUsqr",		
                    word(runTime.endTime().value()),	
                    mesh,			
                    IOobject::MUST_READ,		
                    IOobject::NO_WRITE	
            ),
	    mesh
	); 

	tensor IntegratedMean = Integrate(mean)/timeSpan;
	tensor IntegratedPerb = Integrate(perb)/timeSpan;

	Info << "K(grad U)^2" << endl;
	Info << "-----------" << endl; 
	word component;
	word dir1;
	word dir2;
	
	label c=0;
	for (int i=0 ; i < 3 ; i++ ) {
		switch (i) {
			case 0: 
				dir1 = "x";
				break;
			case 1: 
				dir1 = "y";
				break;
			case 2: 
				dir1 = "z";
				break;
		};
		for (int j=0 ; j < 3 ; j++,c++ ) {
			switch (j) {
				case 0: 
					component="u";
					dir2 = "x";
					break;
				case 1: 
					component="v";
					dir2 = "y";
					break;
				case 2: 
					component="w";
					dir2 = "z";
					break;
			};

			Info << "\t\t d" << component << "/d"<< dir1 << dir1 << " :: " << IntegratedMean.component(c)  << " , " << IntegratedPerb.component(c)  << endl;
					
		} //..for 
	} //.. for
	Info << "-----------------------------------" <<endl << endl;	 
}

void calculateTermBalances::printPressure() { 

	scalar timeSpan = (runTime.endTime() - runTime.startTime()).value();
	getFields("pressure",Scalar)

	scalar IntegratedTotal = Integrate(total)/timeSpan;
	scalar IntegratedMean = Integrate(mean)/timeSpan;
	scalar IntegratedPerb = Integrate(perb)/timeSpan;

	Info << "Potential" << endl;
	Info << "---------" << endl; 
	Info << "total = mean + perb => " << IntegratedTotal  << " = " <<  IntegratedMean << "+ " <<  IntegratedPerb << endl;
	Info << "-----------------------------------" <<endl << endl;	

}

void calculateTermBalances::printNudging() { 

	scalar timeSpan = (runTime.endTime() - runTime.startTime()).value();
	getFields("nudging",Scalar)

	scalar IntegratedTotal = Integrate(total)/timeSpan;
	scalar IntegratedMean = Integrate(mean)/timeSpan;
	scalar IntegratedPerb = Integrate(perb)/timeSpan;

	Info << "Nudging" << endl;
	Info << "---------" << endl; 
	Info << "total = mean + perb => " << IntegratedTotal  << " = " <<  IntegratedMean << "+ " <<  IntegratedPerb << endl;
	Info << "-----------------------------------" <<endl << endl;	

}

void calculateTermBalances::printConvection() { 

	scalar timeSpan = (runTime.endTime() - runTime.startTime()).value();

	getTensor(mean_meanMeanFlux) // mean flux of mean kinetic energy. 
	getTensor(mean_perbPerbFlux) // = Conversion term + reynolds flux along the edge. 
	getTensor(perb_perbMeanFlux) // = Conversion terms (Method II) - only when summed - we do not print it!.
	getTensor(perb_meanPerbFlux) // = Mean Flux of perturbation kinetic energy. 
	getTensor(perb_perbPerbFlux) // = Perturbation flux of perturbation kinetic energy. 

	
	Info << "Mean flux of mean kinetic energy. " << endl;
	Info << "-----------------------------------" <<endl;
	printTensorTerms(mean_meanMeanFlux/timeSpan);
	Info << "-----------------------------------" <<endl << endl;	

	Info << "Conversion term + reynolds flux along the edge. " << endl;
	Info << "-----------------------------------" <<endl;
	printTensorTerms(mean_perbPerbFlux/timeSpan);
	Info << "-----------------------------------" <<endl << endl;	

	Info << "Mean Flux of perturbation kinetic energy.  " << endl;
	Info << "-----------------------------------" <<endl;
	printTensorTerms(perb_meanPerbFlux/timeSpan);
	Info << "-----------------------------------" <<endl << endl;	

	Info << "Perturbation flux of perturbation kinetic energy.  " << endl;
	Info << "-----------------------------------" <<endl;
	printTensorTerms(perb_perbPerbFlux/timeSpan);
	Info << "-----------------------------------" <<endl << endl;	

}

void calculateTermBalances::printConversion() { 

	scalar timeSpan = (runTime.endTime() - runTime.startTime()).value();

	getTensor(mean_perb_conversion)

	Info << "Conversion terms. " << endl;
	Info << "-----------------------------------" <<endl;
	printTensorTerms(mean_perb_conversion/timeSpan);
	Info << "-----------------------------------" <<endl << endl;	
	

}

void calculateTermBalances::printTensorTerms(const tmp<volTensorField>& field) { 
printTensorTerms(field()) ;
}

void calculateTermBalances::printTensorTerms(const volTensorField& field) { 
	word component;
	word dir1;
	word dir2;
	
	label c=0;
	for (int i=0 ; i < 3 ; i++ ) {
		switch (i) {
			case 0: 
				dir1 = "x";
				break;
			case 1: 
				dir1 = "y";
				break;
			case 2: 
				dir1 = "z";
				break;
		};
		for (int j=0 ; j < 3 ; j++,c++ ) {
			switch (j) {
				case 0: 
					component="u";
					dir2 = "x";
					break;
				case 1: 
					component="v";
					dir2 = "y";
					break;
				case 2: 
					component="w";
					dir2 = "z";
					break;
			};

			Info << "\t\t d" << component << "/d"<< dir1 << dir1 << " : " << field.component(c) <<  endl;
					
		} //..for 
	} //.. for 


}


void calculateTermBalances::printAll() { 

	printTemporal(); 
	printConvection();
	printDiffusion();  
	printKUsqr(); 
	printConversion();
	printPotential();
	printNudging();
	printPressure();     


}

