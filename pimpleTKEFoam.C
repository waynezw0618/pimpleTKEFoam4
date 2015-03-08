/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    pimpleFoam

Description
    Large time-step transient solver for incompressible, flow using the PIMPLE
    (merged PISO-SIMPLE) algorithm.

    Sub-models include:
    - turbulence modelling, i.e. laminar, RAS or LES
    - run-time selectable finite volume options, e.g. MRF, explicit porosity

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"
#include "pimpleControl.H"
#include "fvIOoptionList.H"
#include "IOporosityModelList.H"
#include "IOMRFZoneList.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "createFvOptions.H"
    #include "initContinuityErrs.H"

    pimpleControl pimple(mesh);

    //*****************adding reference point of VVCS ************************//
    Info << "searching reference points" <<endl;
    vector refpoint1(6.250420599, 0.02104728876,3.004128853),
           refpoint2(6.250420599, 0.3402290458,3.004128853),
           refpoint3(6.250420599, 0.4335083697,3.004128853),
           refpoint4(6.250420599, 0.7016601976,3.004128853),
           refpoint5(6.250420599, 0.8263088868,3.004128853),     //??? what is
           refpoint6(6.250420599, 0.95414283685,3.004128853),
           refpoint7(6.250420599,0.998026526599,3.004128853),
           refpoint8(6.250420599,0.99950295625,3.004128853);
    
    label refID1(0),refID2(0),refID3(0),refID4(0),refID5(0),refID6(0),refID7(0),refID8(0);
    
    scalar disTol(1e-5);
    
    forAll(U, cellI){
        if(Foam::mag(mesh.C()[cellI]-refpoint1)<disTol){refID1=cellI;}
        if(Foam::mag(mesh.C()[cellI]-refpoint2)<disTol){refID2=cellI;}
        if(Foam::mag(mesh.C()[cellI]-refpoint3)<disTol){refID3=cellI;}
        if(Foam::mag(mesh.C()[cellI]-refpoint4)<disTol){refID4=cellI;}
        if(Foam::mag(mesh.C()[cellI]-refpoint5)<disTol){refID5=cellI;}
        if(Foam::mag(mesh.C()[cellI]-refpoint6)<disTol){refID6=cellI;}
        if(Foam::mag(mesh.C()[cellI]-refpoint7)<disTol){refID7=cellI;}
        if(Foam::mag(mesh.C()[cellI]-refpoint8)<disTol){refID8=cellI;}
    }
    
    Info<<"the 1st refpoint cell["<<refID1<<"] is at "<< mesh.C()[refID1] <<nl
        <<"the 2nd refpoint cell["<<refID2<<"] is at "<< mesh.C()[refID2] <<nl
        <<"the 3rd refpoint cell["<<refID3<<"] is at "<< mesh.C()[refID3] <<nl
        <<"the 4th refpoint cell["<<refID4<<"] is at "<< mesh.C()[refID4] <<nl
        <<"the 5th refpoint cell["<<refID5<<"] is at "<< mesh.C()[refID5] <<nl
        <<"the 6th refpoint cell["<<refID5<<"] is at "<< mesh.C()[refID6] <<nl
        <<"the 7th refpoint cell["<<refID5<<"] is at "<< mesh.C()[refID7] <<nl
        <<"the 8th refpoint cell["<<refID5<<"] is at "<< mesh.C()[refID8] <<endl;
    
    //*********************end reference point *******************************//
    
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

        pU=p*U;
        Uxxx=U.component(vector::X)()*U.component(vector::X)()*U.component(vector::X)();
        Uxxy=U.component(vector::X)()*U.component(vector::X)()*U.component(vector::Y)();
        Uxxz=U.component(vector::X)()*U.component(vector::X)()*U.component(vector::Z)();
        Uyyy=U.component(vector::Y)()*U.component(vector::Y)()*U.component(vector::Y)();
        Uyyx=U.component(vector::Y)()*U.component(vector::Y)()*U.component(vector::X)();
        Uyyz=U.component(vector::Y)()*U.component(vector::Y)()*U.component(vector::Z)();
        Uzzx=U.component(vector::Z)()*U.component(vector::Z)()*U.component(vector::X)();
        Uzzy=U.component(vector::Z)()*U.component(vector::Z)()*U.component(vector::Y)();
        Uzzz=U.component(vector::Z)()*U.component(vector::Z)()*U.component(vector::Z)();

        vorticity=fvc::curl(U);
        Uturb=(U-UMeanMap);
        //ZHANG WEI ADD
        Wturb=(vorticity-vorticityMeanMap);
        const volTensorField gradU(fvc::grad(U));
        ShearStress=0.5*(gradU + gradU.T());
        WpWp=Wturb*Wturb;
      //  Info << "1" <<endl;
        Pens21 =   vorticity.component(vector::X)()*ShearStress.component(tensor::XX)()
                 + vorticity.component(vector::Y)()*ShearStress.component(tensor::YX)()
                 + vorticity.component(vector::Z)()*ShearStress.component(tensor::ZX)();

        Pens22 =   vorticity.component(vector::X)()*ShearStress.component(tensor::XY)()
                 + vorticity.component(vector::Y)()*ShearStress.component(tensor::YY)()
                 + vorticity.component(vector::Z)()*ShearStress.component(tensor::ZY)();

        Pens21 =   vorticity.component(vector::X)()*ShearStress.component(tensor::XZ)()
                 + vorticity.component(vector::Y)()*ShearStress.component(tensor::YZ)()
                 + vorticity.component(vector::Z)()*ShearStress.component(tensor::ZZ)();

        Pens3 = WpWp && ShearStress;

        Pens4 = Uturb*Wturb;

        Tens =  Uturb*magSqr(Wturb);

        DisEnstropy = fvc::grad(Wturb)&&fvc::grad(Wturb)().T();

        // Lambda2

        volTensorField SSplusWW
        (
           (symm(gradU) & symm(gradU)) + (skew(gradU) & skew(gradU))
        );
        //Info << "2" << endl;
        Lambda2=-eigenValues(SSplusWW)().component(vector::Y);


        //epsilon
        epsilon=fvc::grad(Uturb)&&fvc::grad(Uturb)().T();


        Trans= Uturb.component(vector::X)()*Uturb.component(vector::X)()*Uturb.component(vector::Y)()\
               +Uturb.component(vector::Y)()*Uturb.component(vector::Y)()*Uturb.component(vector::Y)()\
               +Uturb.component(vector::Z)()*Uturb.component(vector::Z)()*Uturb.component(vector::Y)();

      epsilonXX=fvc::grad(Uturb)().component(tensor::XX)*fvc::grad(Uturb)().component(tensor::XX) \
                +fvc::grad(Uturb)().component(tensor::XY)*fvc::grad(Uturb)().component(tensor::XY) \
                +fvc::grad(Uturb)().component(tensor::XZ)*fvc::grad(Uturb)().component(tensor::XZ);

      epsilonXY=fvc::grad(Uturb)().component(tensor::XX)*fvc::grad(Uturb)().component(tensor::YX) \
                +fvc::grad(Uturb)().component(tensor::XY)*fvc::grad(Uturb)().component(tensor::YY) \
                +fvc::grad(Uturb)().component(tensor::XZ)*fvc::grad(Uturb)().component(tensor::YZ);

      epsilonXZ=fvc::grad(Uturb)().component(tensor::XX)*fvc::grad(Uturb)().component(tensor::ZX) \
                +fvc::grad(Uturb)().component(tensor::XY)*fvc::grad(Uturb)().component(tensor::ZY) \
                +fvc::grad(Uturb)().component(tensor::XZ)*fvc::grad(Uturb)().component(tensor::ZZ);

      epsilonYY=fvc::grad(Uturb)().component(tensor::YX)*fvc::grad(Uturb)().component(tensor::YX) \
                +fvc::grad(Uturb)().component(tensor::YY)*fvc::grad(Uturb)().component(tensor::YY) \
                +fvc::grad(Uturb)().component(tensor::YZ)*fvc::grad(Uturb)().component(tensor::YZ);

      epsilonYZ=fvc::grad(Uturb)().component(tensor::YX)*fvc::grad(Uturb)().component(tensor::ZX) \
                +fvc::grad(Uturb)().component(tensor::YY)*fvc::grad(Uturb)().component(tensor::ZY) \
                +fvc::grad(Uturb)().component(tensor::YZ)*fvc::grad(Uturb)().component(tensor::ZZ);

      epsilonZZ=fvc::grad(Uturb)().component(tensor::ZX)*fvc::grad(Uturb)().component(tensor::ZX) \
                +fvc::grad(Uturb)().component(tensor::ZY)*fvc::grad(Uturb)().component(tensor::ZY) \
                +fvc::grad(Uturb)().component(tensor::ZZ)*fvc::grad(Uturb)().component(tensor::ZZ);

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
