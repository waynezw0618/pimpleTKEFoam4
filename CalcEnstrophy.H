Info << "Calculating Enstropy ..." <<endl;

const volTensorField gradU(fvc::grad(U));
ShearStress=0.5*(gradU + gradU.T());
vorticity ＝fvc::curl(U)；
Wturb = vorticity - vorticityMeanMap;
WpWp=Wturb*Wturb;
Pens21 = vorticity.component(vector::X)()*ShearStress.component(tensor::XX)()
        + vorticity.component(vector::X)()*ShearStress.component(tensor::XY)()
        + vorticity.component(vector::X)()*ShearStress.component(tensor::XZ)();

Pens22 = vorticity.component(vector::Y)()*ShearStress.component(tensor::YX)()
        + vorticity.component(vector::Y)()*ShearStress.component(tensor::YY)()
        + vorticity.component(vector::Y)()*ShearStress.component(tensor::YZ)();

Pens23 = vorticity.component(vector::Z)()*ShearStress.component(tensor::ZX)()
        + vorticity.component(vector::Z)()*ShearStress.component(tensor::ZY)()
        + vorticity.component(vector::Z)()*ShearStress.component(tensor::ZZ)();

Pens31 = Wturb.component(vector::X)()*Wturb.component(vector::X)()*ShearStress.component(tensor::XX)()
        + Wturb.component(vector::X)()*Wturb.component(vector::Y)()*ShearStress.component(tensor::XY)()
        + Wturb.component(vector::X)()*Wturb.component(vector::Z)()*ShearStress.component(tensor::XZ)();

Pens32 = Wturb.component(vector::Y)()*Wturb.component(vector::X)()*ShearStress.component(tensor::YX)()
        + Wturb.component(vector::Y)()*Wturb.component(vector::Y)()*ShearStress.component(tensor::YY)()
        + Wturb.component(vector::Y)()*Wturb.component(vector::Z)()*ShearStress.component(tensor::YZ)();

Pens33 = Wturb.component(vector::Z)()*Wturb.component(vector::X)()*ShearStress.component(tensor::ZX)()
        + Wturb.component(vector::Z)()*Wturb.component(vector::Y)()*ShearStress.component(tensor::ZY)()
        + Wturb.component(vector::Z)()*Wturb.component(vector::Z)()*ShearStress.component(tensor::ZZ)();


Pens4 = Uturb*Wturb;

Tens1 = Foam::magSqr(Wturb.component(vector::X))*Uturb.component(vector::Y);
Tens2 = Foam::magSqr(Wturb.component(vector::Y))*Uturb.component(vector::Y);
Tens3 = Foam::magSqr(Wturb.component(vector::Z))*Uturb.component(vector::Y);

DisEnstropy1 = Foam::magSqr(fvc::grad(Wturb.component(vector::X)));
DisEnstropy2 = Foam::magSqr(fvc::grad(Wturb.component(vector::Y)));
DisEnstropy3 = Foam::magSqr(fvc::grad(Wturb.component(vector::Z)));

