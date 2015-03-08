Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


volVectorField pU
(
    IOobject
    (
        "pU",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    p*U
);

volScalarField Uxxx
(
    IOobject
    (
        "Uxxx",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    U.component(vector::X)()*U.component(vector::X)()*U.component(vector::X)()
);

volScalarField Uxxy
(
    IOobject
    (
        "Uxxy",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    U.component(vector::X)()*U.component(vector::X)()*U.component(vector::Y)()
);

volScalarField Uxxz
(
    IOobject
    (
        "Uxxz",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    U.component(vector::X)()*U.component(vector::X)()*U.component(vector::Z)()
);

volScalarField Uyyy
(
    IOobject
    (
        "Uyyy",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    U.component(vector::Y)()*U.component(vector::Y)()*U.component(vector::Y)()
);

volScalarField Uyyx
(
    IOobject
    (
        "Uyyx",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    U.component(vector::Y)()*U.component(vector::Y)()*U.component(vector::X)()
);

volScalarField Uyyz
(
    IOobject
    (
        "Uyyz",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    U.component(vector::Y)()*U.component(vector::Y)()*U.component(vector::Z)()
);

volScalarField Uzzz
(
    IOobject
    (
        "Uzzz",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    U.component(vector::Z)()*U.component(vector::Z)()*U.component(vector::Z)()
);

volScalarField Uzzx
(
    IOobject
    (
        "Uzzx",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    U.component(vector::Z)()*U.component(vector::Z)()*U.component(vector::X)()
);

volScalarField Uzzy
(
    IOobject
    (
        "Uzzy",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    U.component(vector::Z)()*U.component(vector::Z)()*U.component(vector::Y)()
);


Info<< "    Calculating vorticity \n" << endl;
volVectorField vorticity
(
   IOobject
   (
        "vorticity",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
   ),
   fvc::curl(U)
);

Info<< "Reading field UMeanMap\n" << endl;
volVectorField UMeanMap
(
    IOobject
    (
        "UMeanMap",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

volVectorField Uturb(U-UMeanMap);

Info<< "Reading field vorticityMeanMap\n" << endl;
volVectorField vorticityMeanMap
(
    IOobject
    (
        "vorticityMeanMap",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

volVectorField Wturb(vorticity-vorticityMeanMap);

Info<< "Reading field URMSMap\n" << endl;
volVectorField URMSMap
(
    IOobject
    (
        "URMSMap",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
     ),
    mesh
);

volTensorField ShearStress
(
    IOobject
    (
        "ShearStress",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    0.5*(fvc::grad(U)+fvc::grad(U)().T())
);

volTensorField WpWp
(
    IOobject
    (
        "ShearStress",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
//    mesh,
//    dimensionedTensor("zero", dimensionSet(0, 0, -2, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0))
    Wturb*Wturb
);

volScalarField Pens21
(
    IOobject
    (
        "Pens21",
         runTime.timeName(),
         mesh,
         IOobject::READ_IF_PRESENT,
         IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(0, 0, -2, 0, 0, 0, 0), 0.0)
);

volScalarField Pens22
(
    IOobject
    (
        "Pens22",
         runTime.timeName(),
         mesh,
         IOobject::READ_IF_PRESENT,
         IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(0, 0, -2, 0, 0, 0, 0), 0.0)
);

volScalarField Pens23
(
    IOobject
    (
        "Pens23",
         runTime.timeName(),
         mesh,
         IOobject::READ_IF_PRESENT,
         IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(0, 0, -2, 0, 0, 0, 0), 0.0)
);

volScalarField Pens3
(
    IOobject
    (
        "Pens3",
         runTime.timeName(),
         mesh,
         IOobject::READ_IF_PRESENT,
         IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(0, 0, -3, 0, 0, 0, 0), 0.0)
);

volTensorField Pens4
(
    IOobject
    (
        "Pens4",
         runTime.timeName(),
         mesh,
         IOobject::READ_IF_PRESENT,
         IOobject::AUTO_WRITE
    ),
//    mesh,
//    dimensionedTensor("zero", dimensionSet(0, 0, -2, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0))
    Uturb*Wturb
);

volVectorField Tens
(
    IOobject
    (
        "Tens",
         runTime.timeName(),
         mesh,
         IOobject::READ_IF_PRESENT,
         IOobject::AUTO_WRITE
    ),
    Uturb*magSqr(Wturb)
 );

volScalarField DisEnstropy
(
    IOobject
    (
        "DisEnstropy",
         runTime.timeName(),
         mesh,
         IOobject::READ_IF_PRESENT,
         IOobject::AUTO_WRITE
    ),
    fvc::grad(Wturb)&&fvc::grad(Wturb)().T()
 );

volScalarField Lambda2
(
    IOobject
    (
        "Lambda2",
         runTime.timeName(),
         mesh,
         IOobject::READ_IF_PRESENT,
         IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(0, 0, -2, 0, 0, 0, 0), 0.0)
);
//epsilon
volScalarField epsilon
(
    IOobject
    (
        "epsilon",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
        fvc::grad(Uturb)&&fvc::grad(Uturb)().T()
);

volScalarField Trans
(
    IOobject
    (
        "Trans",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
        Uturb.component(vector::X)()*Uturb.component(vector::X)()*Uturb.component(vector::Y)()\
       +Uturb.component(vector::Y)()*Uturb.component(vector::Y)()*Uturb.component(vector::Y)()\
       +Uturb.component(vector::Z)()*Uturb.component(vector::Z)()*Uturb.component(vector::Y)()
);

//epsilon
volScalarField epsilonXX
(
    IOobject
    (
        "epsilonXX",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
        fvc::grad(Uturb)().component(tensor::XX)*fvc::grad(Uturb)().component(tensor::XX) \
       +fvc::grad(Uturb)().component(tensor::XY)*fvc::grad(Uturb)().component(tensor::XY) \
       +fvc::grad(Uturb)().component(tensor::XZ)*fvc::grad(Uturb)().component(tensor::XZ)

        //fvc::grad(Uturb.component(vector::X)()).component(vector::X)()
);

volScalarField epsilonXY
(
    IOobject
    (
        "epsilonXY",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
        fvc::grad(Uturb)().component(tensor::XX)*fvc::grad(Uturb)().component(tensor::YX) \
       +fvc::grad(Uturb)().component(tensor::XY)*fvc::grad(Uturb)().component(tensor::YY) \
       +fvc::grad(Uturb)().component(tensor::XZ)*fvc::grad(Uturb)().component(tensor::YZ)
);

volScalarField epsilonXZ
(
    IOobject
    (
        "epsilonXZ",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
        fvc::grad(Uturb)().component(tensor::XX)*fvc::grad(Uturb)().component(tensor::ZX) \
       +fvc::grad(Uturb)().component(tensor::XY)*fvc::grad(Uturb)().component(tensor::ZY) \
       +fvc::grad(Uturb)().component(tensor::XZ)*fvc::grad(Uturb)().component(tensor::ZZ)
);

volScalarField epsilonYY
(
    IOobject
    (
        "epsilonYY",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
        fvc::grad(Uturb)().component(tensor::YX)*fvc::grad(Uturb)().component(tensor::YX) \
       +fvc::grad(Uturb)().component(tensor::YY)*fvc::grad(Uturb)().component(tensor::YY) \
       +fvc::grad(Uturb)().component(tensor::YZ)*fvc::grad(Uturb)().component(tensor::YZ)
);

volScalarField epsilonYZ
(
    IOobject
    (
        "epsilonYZ",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
        fvc::grad(Uturb)().component(tensor::YX)*fvc::grad(Uturb)().component(tensor::ZX) \
       +fvc::grad(Uturb)().component(tensor::YY)*fvc::grad(Uturb)().component(tensor::ZY) \
       +fvc::grad(Uturb)().component(tensor::YZ)*fvc::grad(Uturb)().component(tensor::ZZ)
);

volScalarField epsilonZZ
(
    IOobject
    (
        "epsilonZZ",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
        fvc::grad(Uturb)().component(tensor::ZX)*fvc::grad(Uturb)().component(tensor::ZX) \
       +fvc::grad(Uturb)().component(tensor::ZY)*fvc::grad(Uturb)().component(tensor::ZY) \
       +fvc::grad(Uturb)().component(tensor::ZZ)*fvc::grad(Uturb)().component(tensor::ZZ)
);
//#############VVCS parameters ################//
volScalarField RUlam1
(
    IOobject
    (
        "RUlam1",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
     ),
     Uturb.component(vector::X)()*Lambda2/URMSMap.component(vector::X)
 );

volVectorField RUW1
(
    IOobject
    (
        "RUW1",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    Uturb.component(vector::X)()*Wturb/URMSMap.component(vector::X)
 );

volScalarField RUlam2
(
    IOobject
    (
        "RUlam2",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
     ),
    Uturb.component(vector::X)()*Lambda2/URMSMap.component(vector::X)
 );

volVectorField RUW2
(
    IOobject
    (
        "RUW2",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
     ),
    Uturb.component(vector::X)()*Wturb/URMSMap.component(vector::X)
 );

volScalarField RUlam3
(
    IOobject
    (
        "RUlam3",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
     ),
    Uturb.component(vector::X)()*Lambda2/URMSMap.component(vector::X)
 );

volVectorField RUW3
(
    IOobject
    (
        "RUW3",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
     ),
 Uturb.component(vector::X)()*Wturb/URMSMap.component(vector::X)
 );

volScalarField RUlam4
(
    IOobject
    (
        "RUlam4",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
     ),
    Uturb.component(vector::X)()*Lambda2/URMSMap.component(vector::X)
 );

volVectorField RUW4
(
    IOobject
    (
        "RUW4",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
     ),
    Uturb.component(vector::X)()*Wturb/URMSMap.component(vector::X)
 );

volScalarField RUlam5
(
    IOobject
    (
        "RUlam5",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    Uturb.component(vector::X)()*Lambda2/URMSMap.component(vector::X)
 );

volVectorField RUW5
(
    IOobject
    (
        "RUW5",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
     ),
        Uturb.component(vector::X)()*Wturb/URMSMap.component(vector::X)
 );

volScalarField RUlam6
(
    IOobject
    (
        "RUlam6",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
     ),
    Uturb.component(vector::X)()*Lambda2/URMSMap.component(vector::X)
 );

volVectorField RUW6
(
    IOobject
    (
        "RUW6",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
     ),
    Uturb.component(vector::X)()*Wturb/URMSMap.component(vector::X)
 );

volScalarField RUlam7
(
    IOobject
    (
        "RUlam7",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    Uturb.component(vector::X)()*Lambda2/URMSMap.component(vector::X)
 );

volVectorField RUW7
(
    IOobject
    (
        "RUW7",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
     ),
    Uturb.component(vector::X)()*Wturb/URMSMap.component(vector::X)
 );

volScalarField RUlam8
(
    IOobject
    (
        "RUlam8",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    Uturb.component(vector::X)()*Lambda2/URMSMap.component(vector::X)
 );

volVectorField RUW8
(
    IOobject
    (
        "RUW8",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
     ),
    Uturb.component(vector::X)()*Wturb/URMSMap.component(vector::X)
);
//################END of VVCS ##################/



#include "createPhi.H"


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, mesh.solutionDict().subDict("PIMPLE"), pRefCell, pRefValue);


singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

