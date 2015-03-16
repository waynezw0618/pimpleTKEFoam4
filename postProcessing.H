Info << "post processing for TKE, Enstropy and VVCS" << endl;

Uturb=(U-UMeanMap);

if(TKEbuget)
    #include "CalcTKE.H"

if(Enstrophy)
    #include "CalcEnstrophy.H"

if(VVCS)
    #include "VVCS.H"

/*// Lambda2
volTensorField SSplusWW
(
    (symm(gradU) & symm(gradU)) + (skew(gradU) & skew(gradU))
);

Lambda2=-eigenValues(SSplusWW)().component(vector::Y);
*/