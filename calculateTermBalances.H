#ifndef CALCULATETERMBALANCES_H
#define CALCULATETERMBALANCES_H

#include "dictionary.H"
#include "fvCFD.H"
#include "IDLList.H"

namespace Foam
{


/**
	This function claculates the balances of the various terms after the calculation with the 

*/
class calculateTermBalances {

/****************************************
 *   	Public functions		*
 ****************************************/
public:	
	  calculateTermBalances(fvMesh& imesh,Time& irunTime, word zoneName);

	


/****************************************
 *   	Protected functions		*
 ****************************************/
protected: 

	fvMesh& mesh; 
	Time& runTime;
	volScalarField zoneSelector;
	scalar ZoneVolume;


	// ---------------------------- terms -----------------------------	
	void printTemporal(); 

	void printConvection();
	void printDiffusion();  
	void printDiffusion(); 

	void printPotential();
	void printNudging();
	void printPressure();        

	// ---------------------------- Integration -----------------------------
		scalar Integrate(const volScalarField & field) { 
				scalar sfield = sum(field).value(); 
				reduce(sfield,sumOp<scalar>());
				return sfield;
		}

		scalar Integrate(const tmp<volScalarField> & field) { 
				scalar sfield = sum(field).value(); 
				reduce(sfield,sumOp<scalar>());
				return sfield;
		}

		vector Integrate(const volVectorField & field) { 
				vector sfield = sum(field).value(); 
				reduce(sfield,sumOp<vector>());
				return sfield;
		}


		vector Integrate(const tmp<volVectorField> & field) { 
				vector sfield = sum(field).value(); 
				reduce(sfield,sumOp<vector>());
				return sfield;
		}


		tensor Integrate(const tmp<volTensorField> & field) { 
				tensor sfield = sum(field).value(); 
				reduce(sfield,sumOp<tensor>());
				return sfield;
		}

		tensor Integrate(const volTensorField& field) { 
				tensor sfield = sum(field).value(); 
				reduce(sfield,sumOp<tensor>());
				return sfield;
		}


};


}
#endif 
