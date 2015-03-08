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

        //Q=0.5*(sqr(tr(fvc::grad(U))) - tr(((fvc::grad(U))&(fvc::grad(U)))));
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
