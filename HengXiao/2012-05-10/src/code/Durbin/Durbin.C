#include "Durbin.H"
#include "addToRunTimeSelectionTable.H"
#include "wallFvPatch.H"
#include "fixedInternalValueFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Durbin, 0);
addToRunTimeSelectionTable(RASModel, Durbin, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Durbin::Durbin
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& lamTransportModel
)
:
    RASModel(typeName, U, phi, lamTransportModel),
    GenElliptic(U, phi, lamTransportModel),

    solveK_(coeffDict_.lookupOrAddDefault<Switch>("solveK", true)),
    fBC_(coeffDict_.lookupOrAddDefault<word>("fBC", "automatic")),
    crossTurbDiffusion_(coeffDict_.lookupOrAddDefault<Switch>("crossTurbDiffusion", false)),
    wallsAlignedWithZ_(coeffDict_.lookupOrAddDefault<Switch>("wallsAlignedWithZ", true)),
{
    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Durbin::correct()
{
    GenElliptic::correct();

    if (!turbulence_)
    {
        return;
    }

    volSymmTensorField P = -twoSymm(R_ & fvc::grad(U_));
    volScalarField G("RASModel::G", 0.5*mag(tr(P)));

    volScalarField Ts("T", T());

    #include "../include/epsilonWallI2.H" // set patch internal eps values

    // split R_ into normal diffusion and cross diffusion terms
    volSymmTensorField Rdiag = R_;
    dimensionedScalar kzero = k0_ * 0.0;
    Rdiag.replace(symmTensor::XY, kzero);
    Rdiag.replace(symmTensor::YZ, kzero);
    Rdiag.replace(symmTensor::XZ, kzero);
    volSymmTensorField Rupper = R_ - Rdiag;

    symmTensor  minDiagR = gMin(Rdiag);

    surfaceScalarField Tsf = fvc::interpolate(Ts, "interpolate(T)");
    surfaceSymmTensorField Rdiagf  = fvc::interpolate(Rdiag, "interpolate(R)");
    surfaceSymmTensorField Rupperf = fvc::interpolate(Rupper, "interpolate(R)");

    // Dissipation equation 
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(epsilon_)
        + fvm::div(phi_, epsilon_)
        - fvm::Sp(fvc::div(phi_), epsilon_)
        - fvm::laplacian(Cmu_/sigmaEps_ * Tsf * Rdiagf, epsilon_, "laplacian(epsilon)")
        - fvm::laplacian(nu(), epsilon_, "laplacian(epsilon)")
      ==
        C1_ * G/Ts * ( 1.0 + 0.1*G/epsilon_)
        - fvm::Sp(C2_/Ts, epsilon_)
    );

    if(crossTurbDiffusion_)
    {
        epsEqn() -= fvc::laplacian(Cmu_/sigmaEps_ * Tsf * Rupperf, epsilon_, "laplacian(epsilon)");
    }

    epsEqn().relax();
    epsEqn().boundaryManipulate(epsilon_.boundaryField());
    solve(epsEqn);
    bound(epsilon_, epsilon0_);

    // TKE equation
    if(solveK_)
    {
        tmp<fvScalarMatrix> kEqn
            (
                fvm::ddt(k_)
                + fvm::div(phi_, k_)
                - fvm::Sp(fvc::div(phi_), k_)
                - fvm::laplacian(Cmu_/sigmaK_ * Tsf * Rdiagf, k_, "laplacian(k)")
                - fvm::laplacian(nu(), k_, "laplacian(k)")
                ==
                G
                - fvm::Sp(epsilon_/k_, k_)
            );

        if(crossTurbDiffusion_)
        {
            kEqn() -= fvc::laplacian(Cmu_/sigmaK_ * Tsf * Rupperf, k_, "laplacian(k)");
        }
        
        kEqn().relax();
        solve(kEqn);
    }
    else
    {
        k_ = 0.5 * tr(R_);
    }
    bound(k_, k0_);

    // Reynolds stress equation
    #include "fWallI.H" // set patch internal f values

    tmp<fvSymmTensorMatrix> REqn
        (
            fvm::ddt(R_)
            + fvm::div(phi_, R_)
            - fvm::Sp(fvc::div(phi_), R_)
            - fvm::laplacian(Cmu_/sigmaK_ * Tsf * Rdiagf, R_, "laplacian(R)")
            - fvm::laplacian(nu(), R_, "laplacian(R)")
            + fvm::Sp(epsilon_/k_, R_)
            ==                                        
            P                                         // production tensor
            + k_ * f_
        );

    if(crossTurbDiffusion_)
    {
        REqn() -= fvc::laplacian(Cmu_/sigmaK_*Ts*Rupper, R_, "laplacian(R)");
    }

    REqn().relax();
    solve(REqn);

    if(solveK_)
    {
        forAll(R_, celli)
        {
            symmTensor& rij = R_.internalField()[celli];
            rij.zz() = 2.0*k_.internalField()[celli] - rij.xx() - rij.yy();
        }
    }
    

    volScalarField Ls = L();
    Ts = T(); // re-compute time scale

    volSymmTensorField exSrc = -Clrr1_*dev(R_)/Ts - Clrr2_*dev(P);
    
    
    tmp<fvSymmTensorMatrix> fEqn
        (
            fvm::laplacian(f_)
            ==
            fvm::Sp(1.0/sqr(Ls), f_)
            -
            (
                exSrc/k_ + dev(R_)/(k_*Ts)
            ) / sqr(Ls)
            );
    
    fEqn().relax();
    fEqn().boundaryManipulate(f_.boundaryField());
    solve(fEqn);

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
