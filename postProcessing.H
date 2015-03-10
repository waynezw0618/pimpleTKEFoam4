Info << "post processing for TKE, Enstropy and VVCS" << endl;

Uturb=(U-UMeanMap);
vorticity=fvc::curl(U);
Wturb=(vorticity-vorticityMeanMap);
const volTensorField gradU(fvc::grad(U));
ShearStress=0.5*(gradU + gradU.T());
WpWp=Wturb*Wturb;

＃include "CalcTKE.H"
＃include "CalcEnstrophy.H"

// Lambda2
volTensorField SSplusWW
(
    (symm(gradU) & symm(gradU)) + (skew(gradU) & skew(gradU))
);

Lambda2=-eigenValues(SSplusWW)().component(vector::Y);

#include "VVCS.H"