surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));

volVectorField HbyA("HbyA", U);
HbyA = rAU*UEqn().H();

if (pimple.nCorrPISO() <= 1)
{
    UEqn.clear();
}

//for 2.2.2
surfaceScalarField phiHbyA
(
    "phiHbyA",
    (fvc::interpolate(HbyA) & mesh.Sf())
    + fvc::ddtPhiCorr(rAU, U, phi)
);

//for 2.3.x
/*surfaceScalarField phiHbyA
(
 "phiHbyA",
 (fvc::interpolate(HbyA) & mesh.Sf())
 + rAUf*fvc::ddtCorr(U, phi)
 );*/

fvOptions.relativeFlux(phiHbyA);
//fvOptions.makeRelative(phiHbyA);  //for OF2.3.x

adjustPhi(phiHbyA, U, p);

/* for OF2.3.X
// Update the fixedFluxPressure BCs to ensure flux consistency
setSnGrad<fixedFluxPressureFvPatchScalarField>
(
    p.boundaryField(),
    (
        phiHbyA.boundaryField()
      - fvOptions.relative(mesh.Sf().boundaryField() & U.boundaryField())
    )/(mesh.magSf().boundaryField()*rAUf.boundaryField())
);
*/

// Non-orthogonal pressure corrector loop
while (pimple.correctNonOrthogonal())
{
    // Pressure corrector
    fvScalarMatrix pEqn
    (
        fvm::laplacian(rAUf, p) == fvc::div(phiHbyA)
    );

    pEqn.setReference(pRefCell, pRefValue);

    pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

    if (pimple.finalNonOrthogonalIter())
    {
        phi = phiHbyA - pEqn.flux();
    }
}

#include "continuityErrs.H"

// Explicitly relax pressure for momentum corrector
p.relax();

U = HbyA - rAU*fvc::grad(p);
U.correctBoundaryConditions();
fvOptions.correct(U);
