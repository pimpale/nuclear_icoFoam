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
    Transient solver for incompressible, nuclear laminar flow of Newtonian fluids.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pisoControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[]) {

// clang-format off
#include "setRootCaseLists.H"
#include "createTime.H"
#include "createMesh.H"
// clang-format on

  pisoControl piso(mesh);

  Info << "Reading transportProperties\n" << endl;

  IOdictionary transportProperties(
      IOobject("transportProperties", runTime.constant(), mesh,
               IOobject::MUST_READ_IF_MODIFIED, IOobject::NO_WRITE));

  // Water viscosity
  dimensionedScalar nu("nu", dimViscosity, transportProperties.lookup("nu"));

  // Thermal diffusivity
  dimensionedScalar DT("DT", dimViscosity, transportProperties.lookup("DT"));

  // Diffusivity of thermal neutron concentration
  dimensionedScalar Dnv_t("Dnv_t",
                        dimViscosity,
                        transportProperties.lookup("Dnv_t"));

  // Diffusivity of fast neutron concentration
  dimensionedScalar Dnv_f("Dnv_f",
                        dimViscosity,
                        transportProperties.lookup("Dnv_f"));

  // Neutron multiplier constant for thermal neutrons
  dimensionedScalar nvMul_t("nvMul_t",
                          dimless / dimTime,
                          transportProperties.lookup("nvMul_t"));

  // Neutron multiplier constant for fast neutrons
  dimensionedScalar nvMul_f("nvMul_f",
                          dimless / dimTime,
                          transportProperties.lookup("nvMul_f"));

  // Neutron heating from thermal neutrons
  dimensionedScalar gammaT_t("gammaT_t",
                          dimless / dimTime,
                          transportProperties.lookup("gammaT_t"));

  // Neutron heating from fast neutrons
  dimensionedScalar gammaT_f("gammaT_f",
                          dimless / dimTime,
                          transportProperties.lookup("gammaT_f"));

  // Fast neutron downscattering constant
  dimensionedScalar gamma_f("gamma_f",
                          dimless / dimTime,
                          transportProperties.lookup("gamma_f"));

  // Fast neutron concentration
  Info << "Reading field nv_f\n" << endl;
  volScalarField nv_f(IOobject("nv_f", runTime.timeName(), mesh,
                             IOobject::MUST_READ, IOobject::AUTO_WRITE),
                    mesh);

  // Thermal neutron neutron concentration 
  Info << "Reading field nv_t\n" << endl;
  volScalarField nv_t(IOobject("nv_t", runTime.timeName(), mesh,
                             IOobject::MUST_READ, IOobject::AUTO_WRITE),
                    mesh);

  // Fluid Pressure
  Info << "Reading field p\n" << endl;
  volScalarField p(IOobject("p", runTime.timeName(), mesh, IOobject::MUST_READ,
                            IOobject::AUTO_WRITE),
                   mesh);


  // Fluid velocity
  Info << "Reading field U\n" << endl;
  volVectorField U(IOobject("U", runTime.timeName(), mesh, IOobject::MUST_READ,
                            IOobject::AUTO_WRITE),
                   mesh);

  // Temperature
  Info << "Reading field T\n" << endl;
  volScalarField T(IOobject("T", runTime.timeName(), mesh, IOobject::MUST_READ,
                            IOobject::AUTO_WRITE),
                   mesh);

#include "createPhi.H"

  label pRefCell = 0;
  scalar pRefValue = 0.0;
  setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);
  mesh.setFluxRequired(p.name());
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

    // clang-format off
    fvScalarMatrix nvEqn_t(
          fvm::ddt(nv_t) // over time
        + fvm::div(phi, nv_t) // convection
        - fvm::laplacian(Dnv_t, nv_t) // diffusion
        - nvMul_t * nv_t // fission
        - gamma_f * nv_f // downscattering
    );
    nvEqn_t.solve();

    fvScalarMatrix nvEqn_f(
          fvm::ddt(nv_f)     // over time
        + fvm::div(phi, nv_f) // convection
        - fvm::laplacian(Dnv_f, nv_f) // diffusion
        - nvMul_f * nv_f // fission
        + gamma_f * nv_f // downscattering
      );
    nvEqn_f.solve();

    fvScalarMatrix TEqn(
          fvm::ddt(T) // over time
        + fvm::div(phi, T) // convection
        - fvm::laplacian(DT, T) // diffusion
        + gammaT_t * nv_t // thermal neutron heating
        + gammaT_f * nv_f // fast neutron heating
    );
    TEqn.solve();

    // clang-format on

    runTime.write();

    Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
         << "  ClockTime = " << runTime.elapsedClockTime() << " s" << nl
         << endl;
  }

  Info << "End\n" << endl;

  return 0;
}
