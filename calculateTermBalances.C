#include "calculateTermBalances.H"

#define getFields(name,type)\
	vol##type##Field total(\
            IOobject		\
            (			\
                    "total_energy_#name",		\
                    word(runTime.startTime().value()),	\
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


void calculateTermBalances::printTemporal() { 
	scalar timeSpan = (runTime.endTime() - runTime.startTime()).value();
	getFields("dUdt",Scalar)

	scalar IntegratedTotal = Integrate(total)/timeSpan;
	scalar IntegratedMean = Integrate(mean)/timeSpan;
	scalar IntegratedPerb = Integrate(perb)/timeSpan;

	Info << "Temporal terms" << endl;
	Info << "--------------" << endl;
	Info << "\t total = mean + perb => " << IntegratedTotal << " =  " << IntegratedMean << " + " << IntegratedPerb  << endl;

} //.. printTemporal. 

void calculateTermBalances::printDiffusion() { 
	scalar timeSpan = (runTime.endTime() - runTime.startTime()).value();
	getFields("diffusion",Tensor)

	tensor IntegratedTotal = Integrate(total)/timeSpan;
	tensor IntegratedMean = Integrate(mean)/timeSpan;
	tensor IntegratedPerb = Integrate(perb)/timeSpan;

	word component;
	word dir1;
	word dir2;

	
	Info << "Total Diffusion term (before partial integration)" << endl;
	Info << "-------------------------------------------------" << endl;
	Info << "total = mean + perb" << endl;
	
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

}




