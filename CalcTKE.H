Info<< "Calculating TKE Buget" << endl;

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

//epsilon
epsilon=fvc::grad(Uturb)&&fvc::grad(Uturb)().T();

Trans = Uturb.component(vector::X)()*Uturb.component(vector::X)()*Uturb.component(vector::Y)()\
        +Uturb.component(vector::Y)()*Uturb.component(vector::Y)()*Uturb.component(vector::Y)()\
        +Uturb.component(vector::Z)()*Uturb.component(vector::Z)()*Uturb.component(vector::Y)();

epsilonXX = fvc::grad(Uturb)().component(tensor::XX)*fvc::grad(Uturb)().component(tensor::XX) \
            +fvc::grad(Uturb)().component(tensor::XY)*fvc::grad(Uturb)().component(tensor::XY) \
            +fvc::grad(Uturb)().component(tensor::XZ)*fvc::grad(Uturb)().component(tensor::XZ);

epsilonXY = fvc::grad(Uturb)().component(tensor::XX)*fvc::grad(Uturb)().component(tensor::YX) \
            +fvc::grad(Uturb)().component(tensor::XY)*fvc::grad(Uturb)().component(tensor::YY) \
            +fvc::grad(Uturb)().component(tensor::XZ)*fvc::grad(Uturb)().component(tensor::YZ);

epsilonXZ = fvc::grad(Uturb)().component(tensor::XX)*fvc::grad(Uturb)().component(tensor::ZX) \
            +fvc::grad(Uturb)().component(tensor::XY)*fvc::grad(Uturb)().component(tensor::ZY) \
            +fvc::grad(Uturb)().component(tensor::XZ)*fvc::grad(Uturb)().component(tensor::ZZ);

epsilonYY = fvc::grad(Uturb)().component(tensor::YX)*fvc::grad(Uturb)().component(tensor::YX) \
            +fvc::grad(Uturb)().component(tensor::YY)*fvc::grad(Uturb)().component(tensor::YY) \
            +fvc::grad(Uturb)().component(tensor::YZ)*fvc::grad(Uturb)().component(tensor::YZ);

epsilonYZ = fvc::grad(Uturb)().component(tensor::YX)*fvc::grad(Uturb)().component(tensor::ZX) \
            +fvc::grad(Uturb)().component(tensor::YY)*fvc::grad(Uturb)().component(tensor::ZY) \
            +fvc::grad(Uturb)().component(tensor::YZ)*fvc::grad(Uturb)().component(tensor::ZZ);

epsilonZZ = fvc::grad(Uturb)().component(tensor::ZX)*fvc::grad(Uturb)().component(tensor::ZX) \
            +fvc::grad(Uturb)().component(tensor::ZY)*fvc::grad(Uturb)().component(tensor::ZY) \
            +fvc::grad(Uturb)().component(tensor::ZZ)*fvc::grad(Uturb)().component(tensor::ZZ);
