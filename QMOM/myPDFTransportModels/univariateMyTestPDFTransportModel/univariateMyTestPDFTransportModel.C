/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Code created 2015-2018 by Alberto Passalacqua
    Contributed 2018-07-31 to the OpenFOAM Foundation
    Copyright (C) 2018 OpenFOAM Foundation
    Copyright (C) 2019 Alberto Passalacqua
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

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

#include "univariateMyTestPDFTransportModel.H"
#include "fvcGrad.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PDFTransportModels::univariateMyTestPDFTransportModel
::univariateMyTestPDFTransportModel
(
    const word& name,
    const dictionary& dict,
    const fvMesh& mesh,
    const surfaceScalarField& phi,
    const word& support
)
:
    PDFTransportModel(name, dict, mesh),
    quadrature_(name, mesh, support),
    momentAdvection_
    (
        univariateMomentAdvection::New
        (
            quadrature_.subDict("momentAdvection"),
            quadrature_,
            phi,
            support
        )
    ),
    phim_(phi)
{
    forAll(quadrature_.moments(), momenti)
    {
        mFields_.add(quadrature_.moments()[momenti]);
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::PDFTransportModels::univariateMyTestPDFTransportModel
::~univariateMyTestPDFTransportModel()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::PDFTransportModels::univariateMyTestPDFTransportModel::solve()
{
    //momentAdvection_().update();
    tmp<fv::convectionScheme<scalar>> mvConv
    (
        fv::convectionScheme<scalar>::New
        (
            quadrature_.moments()[0].mesh(), mFields_, phim_, mesh_.divScheme("div(phi.dispersed,moments)")
        )
    );

    // List of moment transport equations
    PtrList<fvScalarMatrix> momentEqns(quadrature_.nMoments());
    
    Info<< "*****Before Solve Equation*****" << endl;
    Info<< "My Test: before Solveing moments equation" << endl;
    forAll (quadrature_.moments()[0], celli)
    {
        if (quadrature_.moments()[0][celli] < -1.0)
        {
            Info<< "m0<1 celli and m0:" << celli << ", " << quadrature_.moments()[0][celli] << endl;
        }    
    }
    
    //********************
    const fvMesh& mesh = this->mesh_;
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    const vectorField& C = mesh.C();
    const scalarField& V = mesh.V();
    
    const surfaceScalarField& CDWeights = mesh.surfaceInterpolation::weights(); 
    
    tmp<volVectorField> tgradM0(fvc::grad(quadrature_.moments()[0]));
    tmp<volVectorField> tgradM1(fvc::grad(quadrature_.moments()[1]));
    tmp<volVectorField> tgradM2(fvc::grad(quadrature_.moments()[2]));
    tmp<volVectorField> tgradM3(fvc::grad(quadrature_.moments()[3]));
    const volVectorField& gradM0 = tgradM0();
    const volVectorField& gradM1 = tgradM1();

    const volVectorField& gradM2 = tgradM2();
    const volVectorField& gradM3 = tgradM3();
    
    tmp<surfaceScalarField> tinterpolationM0(fvc::interpolate(quadrature_.moments()[0]));
    const surfaceScalarField& interpolationM0 = tinterpolationM0();
    
    //forAll(mesh.boundary(), patchi)
    //{
        //Info<< "patch name: " << patchi << endl;
    //}
    
    // output for 202
    /*
    label cells = 82;
    label cellw = 201;
    label cellp = 202;
    label celle = 203;
    label celln = 322;
    */
    // output for 11520
    label cells = 11400;
    label cellw = 11520;
    label cellp = 11521;
    label celle = 11522;
    label celln = 11640;
    
    Info<< "Moments set: " << endl;
    
    Info<< "M0: " //<< endl
    			<< "\'" << cells << "\': " << quadrature_.moments()[0][cells] << ", "
                  	<< "\'" << cellw << "\': " << quadrature_.moments()[0][cellw] << ", "
                 	<< "\'" << cellp << "\': " << quadrature_.moments()[0][cellp] << ", "
                  	<< "\'" << celle << "\': " << quadrature_.moments()[0][celle] << ", "
                  	<< "\'" << celln << "\': " << quadrature_.moments()[0][celln] << endl;
                  
    Info<< "M1: " //<< endl
    			<< "\'" << cells << "\': " << quadrature_.moments()[1][cells] << ", "
                  	<< "\'" << cellw << "\': " << quadrature_.moments()[1][cellw] << ", "
                  	<< "\'" << cellp << "\': " << quadrature_.moments()[1][cellp] << ", "
                  	<< "\'" << celle << "\': " << quadrature_.moments()[1][celle] << ", "
                  	<< "\'" << celln << "\': " << quadrature_.moments()[1][celln] << endl;
                  
    Info<< "M2: " //<< endl
    			<< "\'" << cells << "\': " << quadrature_.moments()[2][cells] << ", "
                  	<< "\'" << cellw << "\': " << quadrature_.moments()[2][cellw] << ", "
                  	<< "\'" << cellp << "\': " << quadrature_.moments()[2][cellp] << ", "
                  	<< "\'" << celle << "\': " << quadrature_.moments()[2][celle] << ", "
                  	<< "\'" << celln << "\': " << quadrature_.moments()[2][celln] << endl;
                  
    Info<< "M3: " //<< endl
    			<< "\'" << cells << "\': " << quadrature_.moments()[3][cells] << ", "
                  	<< "\'" << cellw << "\': " << quadrature_.moments()[3][cellw] << ", "
                  	<< "\'" << cellp << "\': " << quadrature_.moments()[3][cellp] << ", "
                  	<< "\'" << celle << "\': " << quadrature_.moments()[3][celle] << ", "
                  	<< "\'" << celln << "\': " << quadrature_.moments()[3][celln] << endl;
                  
    Info<< "Gradient of moments set: " << endl;
                  
    Info<< "gradM0: " 	<< endl
    			<< "\'" << cells << "\': " << gradM0[cells] << endl
                    	<< "\'" << cellw << "\': " << gradM0[cellw] << endl
                    	<< "\'" << cellp << "\': " << gradM0[cellp] << endl
                    	<< "\'" << celle << "\': " << gradM0[celle] << endl
                    	<< "\'" << celln << "\': " << gradM0[celln] << endl;
                    
    Info<< "gradM1: " 	<< endl
    			<< "\'" << cells << "\': " << gradM1[cells] << endl
                    	<< "\'" << cellw << "\': " << gradM1[cellw] << endl
                    	<< "\'" << cellp << "\': " << gradM1[cellp] << endl
                    	<< "\'" << celle << "\': " << gradM1[celle] << endl
                    	<< "\'" << celln << "\': " << gradM1[celln] << endl;
                    
    Info<< "gradM2: " 	<< endl
    			<< "\'" << cells << "\': " << gradM2[cells] << endl
                    	<< "\'" << cellw << "\': " << gradM2[cellw] << endl
                    	<< "\'" << cellp << "\': " << gradM2[cellp] << endl
                    	<< "\'" << celle << "\': " << gradM2[celle] << endl
                    	<< "\'" << celln << "\': " << gradM2[celln] << endl;
                    
    Info<< "gradM3: " 	<< endl
    			<< "\'" << cells << "\': " << gradM3[cells] << endl
                    	<< "\'" << cellw << "\': " << gradM3[cellw] << endl
                    	<< "\'" << cellp << "\': " << gradM3[cellp] << endl
                    	<< "\'" << celle << "\': " << gradM3[celle] << endl
                    	<< "\'" << celln << "\': " << gradM3[celln] << endl;
                    	
    Info<< "V: " << endl
    			<< "\'" << cells << "\': " << V[cells] << ", "
                  	<< "\'" << cellw << "\': " << V[cellw] << ", "
                  	<< "\'" << cellp << "\': " << V[cellp] << ", "
                  	<< "\'" << celle << "\': " << V[celle] << ", "
                  	<< "\'" << celln << "\': " << V[celln] << endl;

    
    Info<< "owner and neighbour" << endl;
    //label face = 165;
    label face = 22706;
    Info<< " Cell "<< face << ": owner = " << owner[face] << "; neighbour = " << neighbour[face] << endl;
    Info<< "phi[" << face << "] = " << phim_[face] << endl
    	<< "CDWeights[ " << face << "] = " << CDWeights[face] << endl
    	<< "C[" << owner[face] << "]: " << C[owner[face]] << endl 
    	<< "C[" << neighbour[face] << "]: " << C[neighbour[face]] << endl
    	<< "d: " << C[neighbour[face]]  - C[owner[face]] << endl;
    Info<< "interplateM0["<< face << "] = " << interpolationM0[face] << endl;
     
    //face = 401;
    face = 22944;
    Info<< " Cell "<< face << ": owner = " << owner[face] << "; neighbour = " << neighbour[face] << endl;
    Info<< "phi[" << face << "] = " << phim_[face] << endl
    	<< "CDWeights[ " << face << "] = " << CDWeights[face] << endl
    	<< "C[" << owner[face] << "]: " << C[owner[face]] << endl 
    	<< "C[" << neighbour[face] << "]: " << C[neighbour[face]] << endl
    	<< "d: " << C[neighbour[face]]  - C[owner[face]] << endl;
    Info<< "interplateM0["<< face << "] = " << interpolationM0[face] << endl;
    
    //face = 403;
    face = 22945;
    Info<< " Cell "<< face << ": owner = " << owner[face] << "; neighbour = " << neighbour[face] << endl;
    Info<< "phi[" << face << "] = " << phim_[face] << endl
    	<< "CDWeights[ " << face << "] = " << CDWeights[face] << endl
    	<< "C[" << owner[face] << "]: " << C[owner[face]] << endl 
    	<< "C[" << neighbour[face] << "]: " << C[neighbour[face]] << endl
    	<< "d: " << C[neighbour[face]]  - C[owner[face]] << endl;
    Info<< "interplateM0["<< face << "] = " << interpolationM0[face] << endl;
    
    //face = 404;
    face =34486;
    Info<< " Cell "<< face << ": owner = " << owner[face] << "; neighbour = " << neighbour[face] << endl;
    Info<< "phi[" << face << "] = " << phim_[face] << endl
    	<< "CDWeights[ " << face << "] = " << CDWeights[face] << endl
    	<< "C[" << owner[face] << "]: " << C[owner[face]] << endl 
    	<< "C[" << neighbour[face] << "]: " << C[neighbour[face]] << endl
    	<< "d: " << C[neighbour[face]]  - C[owner[face]] << endl;
    Info<< "interplateM0["<< face << "] = " << interpolationM0[face] << endl;
    
    //********************

    // Solve moment transport equations
    forAll(quadrature_.moments(), momenti)
    {
        volScalarMoment& m = quadrature_.moments()[momenti];

        momentEqns.set
        (
            momenti,
            new fvScalarMatrix
            (
                fvm::ddt(m)
              //+ momentAdvection_().divMoments()[momenti]
              //+ fvm::div(phim_, m, "div(phi.dispersed,moments)")
              +  mvConv->fvmDiv(phim_, m)
              ==
                implicitMomentSource(m)
            )
        );
    }
    
    forAll (momentEqns, mEqni)
    {
        //Info<< momentEqns[mEqni] << endl;
        momentEqns[mEqni].relax();
        momentEqns[mEqni].solve();
    }


    Info<< "*****After Solve Equation*****" << endl;
    Info<< "My Test: before updateQuadrature" << endl;
    forAll (quadrature_.moments()[0], celli)
    {
        if (quadrature_.moments()[0][celli] > 1.0)
        {
            if (quadrature_.moments()[1][celli] < 0.0)
            {
                Info<< "m0>1 and m1<0, celli:" << celli << endl;
                Info<< "M0 and M1 and M2: " << quadrature_.moments()[0][celli] 
                    << ", " << quadrature_.moments()[1][celli]
                    << ", " << quadrature_.moments()[2][celli]
                    << endl;
            }
        }
        if (quadrature_.moments()[0][celli] < -1.0)
        {
            Info<< "m0<1 celli and m0:" << celli << ", " << quadrature_.moments()[0][celli] << endl;
        }
        
    }
    
    Info<< "Moments set: " << endl;
    
    Info<< "M0: " //<< endl
    			<< "\'" << cells << "\': " << quadrature_.moments()[0][cells] << ", "
                  	<< "\'" << cellw << "\': " << quadrature_.moments()[0][cellw] << ", "
                 	<< "\'" << cellp << "\': " << quadrature_.moments()[0][cellp] << ", "
                  	<< "\'" << celle << "\': " << quadrature_.moments()[0][celle] << ", "
                  	<< "\'" << celln << "\': " << quadrature_.moments()[0][celln] << endl;
                  
    Info<< "M1: " //<< endl
    			<< "\'" << cells << "\': " << quadrature_.moments()[1][cells] << ", "
                  	<< "\'" << cellw << "\': " << quadrature_.moments()[1][cellw] << ", "
                  	<< "\'" << cellp << "\': " << quadrature_.moments()[1][cellp] << ", "
                  	<< "\'" << celle << "\': " << quadrature_.moments()[1][celle] << ", "
                  	<< "\'" << celln << "\': " << quadrature_.moments()[1][celln] << endl;
                  
    Info<< "M2: " //<< endl
    			<< "\'" << cells << "\': " << quadrature_.moments()[2][cells] << ", "
                  	<< "\'" << cellw << "\': " << quadrature_.moments()[2][cellw] << ", "
                  	<< "\'" << cellp << "\': " << quadrature_.moments()[2][cellp] << ", "
                  	<< "\'" << celle << "\': " << quadrature_.moments()[2][celle] << ", "
                  	<< "\'" << celln << "\': " << quadrature_.moments()[2][celln] << endl;
                  
    Info<< "M3: " //<< endl
    			<< "\'" << cells << "\': " << quadrature_.moments()[3][cells] << ", "
                  	<< "\'" << cellw << "\': " << quadrature_.moments()[3][cellw] << ", "
                  	<< "\'" << cellp << "\': " << quadrature_.moments()[3][cellp] << ", "
                  	<< "\'" << celle << "\': " << quadrature_.moments()[3][celle] << ", "
                  	<< "\'" << celln << "\': " << quadrature_.moments()[3][celln] << endl;

    
    // output for 183  
/*
    
    Info<< "Moments set: " << endl;
    Info<< "M0[celli=63]: " << quadrature_.moments()[0][63] << endl;
    Info<< "M0[celli=182]: " << quadrature_.moments()[0][182] << endl;
    Info<< "M0[celli=183]: " << quadrature_.moments()[0][183] << endl;
    Info<< "M0[celli=184]: " << quadrature_.moments()[0][184] << endl;
    Info<< "M0[celli=303]: " << quadrature_.moments()[0][303] << endl;
    
    Info<< "M1[celli=63]: " << quadrature_.moments()[1][63] << endl;
    Info<< "M1[celli=182]: " << quadrature_.moments()[1][182] << endl;
    Info<< "M1[celli=183]: " << quadrature_.moments()[1][183] << endl;
    Info<< "M1[celli=184]: " << quadrature_.moments()[1][184] << endl;
    Info<< "M1[celli=303]: " << quadrature_.moments()[1][303] << endl;
    
    Info<< "M2[celli=63]: " << quadrature_.moments()[2][63] << endl;
    Info<< "M2[celli=182]: " << quadrature_.moments()[2][182] << endl;
    Info<< "M2[celli=183]: " << quadrature_.moments()[2][183] << endl;
    Info<< "M2[celli=184]: " << quadrature_.moments()[2][184] << endl;
    Info<< "M2[celli=303]: " << quadrature_.moments()[2][303] << endl;
    
    Info<< "M3[celli=63]: " << quadrature_.moments()[3][63] << endl;
    Info<< "M3[celli=182]: " << quadrature_.moments()[3][182] << endl;
    Info<< "M3[celli=183]: " << quadrature_.moments()[3][183] << endl;
    Info<< "M3[celli=184]: " << quadrature_.moments()[3][184] << endl;
    Info<< "M3[celli=303]: " << quadrature_.moments()[3][303] << endl;
    
    Info<< "Phi: " << endl;
    Info<< "phi[facei=127]: " << phim_[127] << endl;
    Info<< "phi[facei=363]: " << phim_[363] << endl;
    Info<< "phi[facei=365]: " << phim_[365] << endl;
    Info<< "phi[facei=366]: " << phim_[366] << endl;
    
    Info<< "owner and neighbour" << endl;
    label face = 127;
    Info<< " Cell "<< face << ": owner = " << owner[face] << "; neighbour = " << neighbour[face] << endl
    << "CDWeights[" << face << "] = " << CDWeights[face] << endl
    << "gradM0: own = " << gradM0[owner[face]] << "; neighbour = " << gradM0[neighbour[face]] << endl
    << "gradM1: own = " << gradM1[owner[face]] << "; neighbour = " << gradM1[neighbour[face]] << endl
    << "gradM2: own = " << gradM2[owner[face]] << "; neighbour = " << gradM2[neighbour[face]] << endl
    << "gradM3: own = " << gradM3[owner[face]] << "; neighbour = " << gradM3[neighbour[face]] << endl
    << "C[" << owner[face] << "]: " << C[owner[face]] << "; C[" << neighbour[face] << "]: " << C[neighbour[face]] << endl
    << "d: " << C[neighbour[face]]  - C[owner[face]] << endl;  
    
    Info<< "interplateM0["<< face << "] = " << interpolationM0[face] << endl;
     
    face = 363;
    Info<< " Cell "<< face << ": owner = " << owner[face] << "; neighbour = " << neighbour[face] << endl
    << "CDWeights[ " << face << "] = " << CDWeights[face] << endl
    << "gradM0: own = " << gradM0[owner[face]] << "; neighbour = " << gradM0[neighbour[face]] << endl
    << "gradM1: own = " << gradM1[owner[face]] << "; neighbour = " << gradM1[neighbour[face]] << endl
    << "gradM2: own = " << gradM2[owner[face]] << "; neighbour = " << gradM2[neighbour[face]] << endl
    << "gradM3: own = " << gradM3[owner[face]] << "; neighbour = " << gradM3[neighbour[face]] << endl
    << "C[" << owner[face] << "]: " << C[owner[face]] << "; C[" << neighbour[face] << "]: " << C[neighbour[face]] << endl
    << "d: " << C[neighbour[face]]  - C[owner[face]] << endl;
    
    Info<< "interplateM0["<< face << "] = " << interpolationM0[face] << endl;
    
    face = 365;
    Info<< " Cell "<< face << ": owner = " << owner[face] << "; neighbour = " << neighbour[face] << endl
    << "CDWeights[ " << face << "] = " << CDWeights[face] << endl
    << "gradM0: own = " << gradM0[owner[face]] << "; neighbour = " << gradM0[neighbour[face]] << endl
    << "gradM1: own = " << gradM1[owner[face]] << "; neighbour = " << gradM1[neighbour[face]] << endl
    << "gradM2: own = " << gradM2[owner[face]] << "; neighbour = " << gradM2[neighbour[face]] << endl
    << "gradM3: own = " << gradM3[owner[face]] << "; neighbour = " << gradM3[neighbour[face]] << endl
    << "C[" << owner[face] << "]: " << C[owner[face]] << "; C[" << neighbour[face] << "]: " << C[neighbour[face]] << endl
    << "d: " << C[neighbour[face]]  - C[owner[face]] << endl;
    
    Info<< "interplateM0["<< face << "] = " << interpolationM0[face] << endl;
    
    face = 366;
    Info<< " Cell "<< face << ": owner = " << owner[face] << "; neighbour = " << neighbour[face] << endl
    << "CDWeights[ " << face << "] = " << CDWeights[face] << endl
    << "gradM0: own = " << gradM0[owner[face]] << "; neighbour = " << gradM0[neighbour[face]] << endl
    << "gradM1: own = " << gradM1[owner[face]] << "; neighbour = " << gradM1[neighbour[face]] << endl
    << "gradM2: own = " << gradM2[owner[face]] << "; neighbour = " << gradM2[neighbour[face]] << endl
    << "gradM3: own = " << gradM3[owner[face]] << "; neighbour = " << gradM3[neighbour[face]] << endl
    << "C[" << owner[face] << "]: " << C[owner[face]] << "; C[" << neighbour[face] << "]: " << C[neighbour[face]] << endl
    << "d: " << C[neighbour[face]]  - C[owner[face]] << endl;
    
    Info<< "interplateM0["<< face << "] = " << interpolationM0[face] << endl;
    
*/
/*
    Info<< "M0: " << quadrature_.moments()[0] << endl;

    Info<< "M0[celli=20894]: " << quadrature_.moments()[0][20894] << endl;
    Info<< "M0[celli=21093]: " << quadrature_.moments()[0][21093] << endl;
    Info<< "M0[celli=21094]: " << quadrature_.moments()[0][21094] << endl;
    Info<< "M0[celli=21095]: " << quadrature_.moments()[0][21095] << endl;
    Info<< "M0[celli=21294]: " << quadrature_.moments()[0][21294] << endl;
    Info<< "phi[celli=41685]: " << phim_[41685] << endl;
    Info<< "phi[celli=42081]: " << phim_[42081] << endl;
    Info<< "phi[celli=42083]: " << phim_[42083] << endl;
    Info<< "phi[celli=42084]: " << phim_[42084] << endl;
*/    
    quadrature_.updateQuadrature();

    Info<< "My Test: after updateQuadrature" << endl;
    //Info<< "Perform univariatePDF: " << endl;

    evaporation();

    if (solveMomentSources())
    {
        this->explicitMomentSource();
    }
}


// ************************************************************************* //
