Info << "Calculating Enstropy ..." <<endl;
Pens21 = vorticity.component(vector::X)()*ShearStress.component(tensor::XX)()
        + vorticity.component(vector::Y)()*ShearStress.component(tensor::YX)()
        + vorticity.component(vector::Z)()*ShearStress.component(tensor::ZX)();

Pens22 = vorticity.component(vector::X)()*ShearStress.component(tensor::XY)()
        + vorticity.component(vector::Y)()*ShearStress.component(tensor::YY)()
        + vorticity.component(vector::Z)()*ShearStress.component(tensor::ZY)();

Pens21 = vorticity.component(vector::X)()*ShearStress.component(tensor::XZ)()
        + vorticity.component(vector::Y)()*ShearStress.component(tensor::YZ)()
        + vorticity.component(vector::Z)()*ShearStress.component(tensor::ZZ)();

Pens3 = WpWp && ShearStress;

Pens4 = Uturb*Wturb;

Tens = Uturb*magSqr(Wturb);

DisEnstropy = fvc::grad(Wturb)&&fvc::grad(Wturb)().T();