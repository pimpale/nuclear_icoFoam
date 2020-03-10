/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    icoFoam

Description
    Transient solver for incompressible, laminar flow of Newtonian fluids.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pisoControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[]) {
#include "createMesh.H"
#include "createTime.H"
#include "setRootCaseLists.H"

  pisoControl piso(mesh);

  // BEGIN createFields

  Info << "Reading transportProperties\n" << endl;

  IOdictionary transportProperties(
      IOobject("transportProperties", runTime.constant(), mesh,
               IOobject::MUST_READ_IF_MODIFIED, IOobject::NO_WRITE));

  dimensionedScalar nu("nu", dimViscosity, transportProperties.lookup("nu"));

  // dimensionedScalar

  // Macroscopic nuclear absorption cross section
  dimensionedScalar sigmaA("sigmaA",
                           dimless / dimLength, // (m -1)
                           transportProperties.lookup("sigmaA"));

  // Macroscopic nuclear fission cross section
  dimensionedScalar sigmaF("sigmaF",
                           dimless / dimLength, // (m -1)
                           transportProperties.lookup("sigmaF"));


  // Macroscopic nuclear scatter cross section
  dimensionedScalar sigmaS("sigmaS",
                           dimless / dimLength, // (m -1)
                           transportProperties.lookup("sigmaS"));

  Info << "Reading field p\n" << endl;
  volScalarField p(IOobject("p", runTime.timeName(), mesh, IOobject::MUST_READ,
                            IOobject::AUTO_WRITE),
                   mesh);

  Info << "Reading field nPhi\n" << endl;
  volScalarField nv(IOobject("nPhi", runTime.timeName(), mesh,
                             IOobject::MUST_READ, IOobject::AUTO_WRITE),
                    mesh);

  Info << "Reading field U\n" << endl;
  volVectorField U(IOobject("U", runTime.timeName(), mesh, IOobject::MUST_READ,
                            IOobject::AUTO_WRITE),
                   mesh);

#include "createPhi.H"

  label pRefCell = 0;
  scalar pRefValue = 0.0;
  setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);
  mesh.setFluxRequired(p.name());
  
  // END createFields

#include "initContinuityErrs.H"

  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

  Info << "\nStarting time loop\n" << endl;

  while (runTime.loop()) {
    Info << "Time = " << runTime.timeName() << nl << endl;

#include "CourantNo.H"

    // Momentum predictor

    fvVectorMatrix UEqn(fvm::ddt(U) + fvm::div(phi, U) - fvm::laplacian(nu, U));

    if (piso.momentumPredictor()) {
      solve(UEqn == -fvc::grad(p));
    }

    // --- PISO loop
    while (piso.correct()) {
      volScalarField rAU(1.0 / UEqn.A());
      volVectorField HbyA(constrainHbyA(rAU * UEqn.H(), U, p));
      surfaceScalarField phiHbyA("phiHbyA",
                                 fvc::flux(HbyA) + fvc::interpolate(rAU) *
                                                       fvc::ddtCorr(U, phi));

      adjustPhi(phiHbyA, U, p);

      // Update the pressure BCs to ensure flux consistency
      constrainPressure(p, U, phiHbyA, rAU);

      // Non-orthogonal pressure corrector loop
      while (piso.correctNonOrthogonal()) {
        // Pressure corrector

        fvScalarMatrix pEqn(fvm::laplacian(rAU, p) == fvc::div(phiHbyA));

        pEqn.setReference(pRefCell, pRefValue);

        pEqn.solve();

        if (piso.finalNonOrthogonalIter()) {
          phi = phiHbyA - pEqn.flux();
        }
      }

#include "continuityErrs.H"

      U = HbyA - rAU * fvc::grad(p);
      U.correctBoundaryConditions();
    }

    // add these lines...
    fvScalarMatrix nvEqn(fvm::ddt(nv) + fvm::div(phi, nv) -
                         fvm::laplacian(Dnv, nv) - nvMul * nv);

    nvEqn.solve();
    // done adding lines...

    runTime.write();

    Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
         << "  ClockTime = " << runTime.elapsedClockTime() << " s" << nl
         << endl;
  }

  Info << "End\n" << endl;

  return 0;
}

// ************************************************************************* //
