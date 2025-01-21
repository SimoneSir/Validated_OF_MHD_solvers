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

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    pisoControl piso(mesh);

    #include "createFields.H"
    #include "initContinuityErrs.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "CourantNo.H"

		for (int corr=0; corr<PotEcorr; corr++)	//	MHD coupling loop
		{
            // Momentum predictor

            fvVectorMatrix UEqn
            (
                fvm::ddt(U)
            + fvm::div(phi, U)
            - fvm::laplacian(nu, U)
            ==
            (1.0/rho) * lorentz		// MHD part
            );

            if (piso.momentumPredictor())
            {
                solve(UEqn == -fvc::grad(p));
            }

            // --- PISO loop
            while (piso.correct())
            {
                volScalarField rAU(1.0/UEqn.A());
                volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
                surfaceScalarField phiHbyA
                (
                    "phiHbyA",
                    fvc::flux(HbyA)
                + fvc::interpolate(rAU)*fvc::ddtCorr(U, phi)
                );

                adjustPhi(phiHbyA, U, p);

                // Update the pressure BCs to ensure flux consistency
                constrainPressure(p, U, phiHbyA, rAU);

                // Non-orthogonal pressure corrector loop
                while (piso.correctNonOrthogonal())
                {
                    // Pressure corrector

                    fvScalarMatrix pEqn
                    (
                        fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
                    );

                    pEqn.setReference(pRefCell, pRefValue);
                    pEqn.solve();

                    if (piso.finalNonOrthogonalIter())
                    {
                        phi = phiHbyA - pEqn.flux();
                    }
                }

                #include "continuityErrs.H"

                U = HbyA - rAU*fvc::grad(p);
                U.correctBoundaryConditions();
            }


            // MHD part
            //surfaceScalarField xi = sigma * fvc::interpolate(U ^ B) & mesh.Sf();
            surfaceScalarField xi = sigma * (fvc::interpolate(U) ^ B) & mesh.Sf(); 
   		
            fvScalarMatrix PotEEqn
            (
                fvm::laplacian(sigma,PotE) == fvc::div(xi)
            );

            PotEEqn.setReference(PotERefCell, PotERefValue);
            PotEEqn.solve();
            //PotE.correctBoundaryConditions();
		
		    surfaceScalarField jn = -sigma * fvc::snGrad(PotE) * mesh.magSf() + xi;
            surfaceVectorField jnv = jn * mesh.Cf();
		
		    jcenter = fvc::surfaceIntegrate(jnv) - fvc::surfaceIntegrate(jn) * mesh.C();
            //jcenter.correctBoundaryConditions();

		    //lorentz = jcenter ^ B;
            lorentz = -fvc::surfaceIntegrate(jn * ( B ^ mesh.Cf() ) ) - (mesh.C() ^ fvc::surfaceIntegrate(jn * B));
            //lorentz.correctBoundaryConditions();

        }
        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}
// ************************************************************************* //