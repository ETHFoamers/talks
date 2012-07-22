#include "muGenElliptic.H"
#include "wallFvPatch.H"
#include "wallDistData.H"
#include "wallPointYPlus.H"
#include "gaussLaplacianScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
  makeFvLaplacianTypeScheme(gaussLaplacianScheme, symmTensor, symmTensor)
}
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

muGenElliptic::muGenElliptic
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& lamTransportModel,
    const volSymmTensorField & RAvg,
    const volScalarField & epsilonAvg
)
:
    muRASModel(typeName, U, phi, lamTransportModel, RAvg, epsilonAvg),
    mesh_(U.mesh()),
  
{
    IOdictionary relaxParameters
        (
            IOobject
            (
                "relaxParameters",
                runTime_.constant(),
                "../../constant",
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

    dictionary couplingDict(relaxParameters.subDictPtr("couplingOptions"));

    imposeTurbEvery_ = 
        couplingDict.lookupOrDefault<label>("mapL2REvery", 1, true);
    
    Info << "(If enabled) directly imposing turb quantities every: " << imposeTurbEvery_ << endl;
}


void muGenElliptic::updateKolmogorovFlag()
{
    // wall unit as defined by nu/sqrt(tauw/rho)
    volScalarField ystar
    (
        IOobject
        (
            "ystar",
            mesh_.time().constant(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("ystar", dimLength, 1.0)
    );

    const fvPatchList& patches = mesh_.boundary();
    forAll(patches, patchi)
    {
        if (isA<wallFvPatch>(patches[patchi]))
        {
            const fvPatchVectorField& Uw = U_.boundaryField()[patchi];
            const scalarField& nuw = nu().boundaryField()[patchi];
            // Note: nuw is used instead of nueff
            // for wall-resolving mesh, nut should be zero at wall
            ystar.boundaryField()[patchi] =
                nuw/sqrt(nuw*mag(Uw.snGrad()) + VSMALL);
        }
    }

    wallPointYPlus::yPlusCutOff = 500;
    wallDistData<wallPointYPlus> y(mesh_, ystar);
    
    KolmogorovFlag_ = pos(yStarLim_ - y/ystar);
}


tmp<volScalarField> muGenElliptic::T() const
{
    return max
        (
            k_/(epsilon_ + epsilonSmall_),
            KolmogorovFlag_ * 6.0 * sqrt(nu()/(epsilon_ + epsilonSmall_))
        );
}

tmp<volScalarField> muGenElliptic::L() const
{
    return
        CL_*max
        (
            pow(k_,1.5)/(epsilon_ + epsilonSmall_),
            KolmogorovFlag_ * CEta_ * pow(pow(nu(),3.0)/(epsilon_ + epsilonSmall_),0.25)
        );
}


tmp<volSymmTensorField> muGenElliptic::devReff() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "devRhoReff",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            R_ - nu()*dev(twoSymm(fvc::grad(U_)))
        )
    );
}


tmp<fvVectorMatrix> muGenElliptic::divDevReff(volVectorField& U) const
{
    if(implicitDiv_)
    {
        return
            (
                fvc::div(R_)
                + fvc::laplacian(nut(), U, "laplacian(nuEff,U)")
                - fvm::laplacian(nuEff(), U)
            );
    }
    else
    {
        return
            (
                fvc::div(R_)
                - fvm::laplacian(nu(), U)
            );
    }
}


void muGenElliptic::correct()
{
    updateKolmogorovFlag();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
