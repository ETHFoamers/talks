#include "mcParticleCloud.H"
#include "fvMesh.H"
#include "interpolationCellPointFace.H"
#include "boundBox.H"
#include "fvc.H"

namespace Foam
{

    defineTemplateTypeNameAndDebug(Cloud<mcParticle>, 0);
}


Foam::mcParticleCloud::mcParticleCloud
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& cloudName,
    const volVectorField* U,
    const volScalarField* p,
    volScalarField* rho
)
:
    Cloud<mcParticle>(mesh, cloudName, false),
    mesh_(mesh),
    dict_(dict),
    runTime_(mesh.time()),
    Ufv_
    (
        U ? *U : mesh_.lookupObject<volVectorField>
                     (dict_.lookupOrDefault<word>("UName", "U"))
    ),
    pfv_
    (
        p ? *p : mesh_.lookupObject<volScalarField>
                     (dict_.lookupOrDefault<word>("pName", "p"))
    ),
    mMom_
    (
        IOobject
        (
            "mMoment",
            runTime_.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("mMoment", dimMass, 0)
    ),
    VMom_
    (
        IOobject
        (
            "VMoment",
            runTime_.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("VMoment", dimVolume, 0)
    ),
    UMom_
    (
        IOobject
        (
            "UMoment",
            runTime_.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("UMoment", dimMass*dimVelocity, vector::zero)
    ),
    uuMom_
    (
        IOobject
        (
            "uuMoment",
            runTime_.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("uuMoment", dimEnergy, symmTensor::zero)
    ),
    rho_
    (
        rho ? *rho : const_cast<volScalarField&>(
            mesh_.lookupObject<volScalarField>(
                dict_.lookupOrDefault<word>("rhoName", "rho")))
    ),
    pmd_
    (
        IOobject
        (
            "pmd",
            runTime_.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimDensity,
        mMom_/mesh_.V(),
        rhocPdf_.boundaryField()
    ),
    U_
    (
        IOobject
        (
            "UPdf",
            runTime_.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimVelocity,
        UMom_/max(mMom_, SMALL_MASS),
        Ufv_.boundaryField()
    ),
    Tau_
    (
        IOobject
        (
            "TauPdf",
            runTime_.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        /* SKIPPED INIT AND BC */
    ),

    k_
    (
        IOobject
        (
            "kPdf",
            runTime_.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimVelocity*dimVelocity,
        0.5*tr(TaucPdf_.dimensionedInternalField()),
        /* SKIPPED BC */
    )
{
    // If particle data found, read from files. Otherwise, initialize.
    if (size() > 0)
    {
        mcParticle::readFields(*this);
    }
    else
    {
        initReleaseParticles();
    }
    initMoments();
    updateParticlePDF();
}


Foam::scalar Foam::mcParticleCloud::evolve()
{
    // SKIPPED RELEASE PARTICLES AT INLET PATCHES

    // Time-averaging factor
    scalar existWt = 1.0/(1.0 + (runTime_.deltaT()/AvgTimeScale_).value());

    mcParticle::trackingData td(/* ... */);

    Cloud<mcParticle>::move(td, runTime_.deltaT().value());

    // Extract statistical averaging to obtain mesh-based quantities
    updateCloudPDF(existWt);
}


void Foam::mcParticleCloud::updateCloudPDF(scalar existWt)
{
    DimensionedField<scalar, volMesh> mMomInstant
    (
        IOobject
        (
            "mMomInstant",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("mMomInstant", dimMass, 0.0)
    );
    // ...
    // Similar for VMomInstant, UMomInstant and uuMomInstant
    // ...

    // Loop through particles to accumulate moments (0, 1, 2 order)
    interpolationCellPointFace<vector> UInterp(Ufv_);
    forAllConstIter(mcParticleCloud, *this, pIter)
    {
        const mcParticle& p = pIter();
        label cellI = p.cell();
        vector U = UInterp.interpolate(p.position(), cellI, p.face());
        vector u = p.U() - U;

        const scalar m = p.m();
        mMomInstant[cellI]  += m;
        VMomInstant[cellI]  += m / p.rho();
        UMomInstant[cellI]  += m * p.U();
        uuMomInstant[cellI] += m * symm(u*u);
    }

    scalar newWt = 1.0 - existWt;
    // Do time-averaging of moments and compute mean fields
    mMom_  = existWt * mMom_  + newWt * mMomInstant;
    pmd_.internalField() = mMom_ / mesh_.V();

    VMom_  = existWt * VMom_  + newWt * VMomInstant;
    rho_.internalField()   = mMom_ / VMom_;
    rho_.correctBoundaryConditions();

    UMom_  = existWt * UMom_  + newWt * UMomInstant;
    U_.internalField()   = UMom_ / mMom_;
    U_.correctBoundaryConditions();

    uuMom_ = existWt * uuMom_ + newWt * uuMomInstant;
    Tau_.internalField() = uuMom_/mMom_;
    Tau_.correctBoundaryConditions();

    k_.internalField()   = 0.5*tr(Tau_.internalField());
    k_.correctBoundaryConditions();
}


void Foam::mcParticleCloud::initReleaseParticles()
{
    // Populate each cell with 30 particles in each cell
    forAll(Ufv_, celli)
    {
        scalar m = mesh_.V()[celli]*rho_[celli] / 30;
        vector U = U_[celli];
        scalar urms = sqrt(2./3.*kfv()()[celli]);
        vector uscales(urms, urms, urms);
        particleGenInCell(celli, 30, m, U, uscales);
    }
}


void Foam::mcParticleCloud::particleGenInCell
(
    label celli,
    label N,
    scalar m,
    const vector& Updf,
    const vector& uscales
)
{
    boundBox cellbb
    (
        pointField
        (
            mesh_.points(),
            mesh_.cellPoints()[celli]
        ),
        false
    );

    vector minb = cellbb.min();
    vector dimb = cellbb.max() - minb;

    label Npgen = 0;
    for (int i = 0; i < 100*N; ++i)
    {
        // Relative coordinate [0, 1] in this cell
        vector xi = random().vector01();
        // Random offset from min point
        scalar rx = min(max(10.0*SMALL, xi.x()), 1.0-10.0*SMALL);
        scalar ry = min(max(10.0*SMALL, xi.y()), 1.0-10.0*SMALL);
        scalar rz = min(max(10.0*SMALL, xi.z()), 1.0-10.0*SMALL);
        vector offsetRnd(rx*dimb.x(), ry*dimb.y(), rz*dimb.z());

        // Generate a particle position
        vector position = minb + offsetRnd;

        // If the case has reduced dimensionality, put the coordinate of the
        // reduced dimension onto the coordinate plane
        if (mesh_.nGeometricD() <= 2)
        {
            meshTools::constrainDirection(mesh_, mesh_.geometricD(), position);
        }

        // Initially put N particle per cell
        if (mesh_.pointInCell(position, celli))
        {
            // random() not shown here
            vector u
            (
              random().GaussNormal()*uscales.x(),
              random().GaussNormal()*uscales.y(),
              random().GaussNormal()*uscales.z()
            );
            vector Up = u + U;

            mcParticle* ptr = new mcParticle
            (
                *this,
                position,
                celli,
                m,
                Up,
                rho_[celli]
            );

            addParticle(ptr);
            ++Npgen;
        }

        // until enough particles are generated.
        if (Npgen >= N) break;
    }

    if (Npgen < N)
    {
        FatalErrorIn("mcParticleCloud::initReleaseParticles()")
            << "Only " << Npgen << " particles generated for cell "
            << celli << nl << "Something is wrong" << exit(FatalError);
    }
}


void Foam::mcParticleCloud::initMoments()
{
    bool readOk =
        mMom_.headerOk() &&
        VMom_.headerOk() &&
        UMom_.headerOk() &&
        uuMom_.headerOk();
    if (readOk)
    {
        Info<< "Moments read correctly." << endl;
    }
    else if (size() > 0)
    {
        Info<< "Moments are missing. Forced re-initialization." << endl;
        mMom_  = mesh_.V()*rho_;
        pmd_.internalField() = rho_.internalField();

        VMom_  = mMom_/rho_;

        UMom_  = mMom_*Ufv_;
        U_.internalField() = Ufv_.internalField();
        U_.correctBoundaryConditions();

        uuMom_ = mMom_*turbulenceModel().R()();

        k_.internalField() = kfv().internalField();
        k_.correctBoundaryConditions();
    }
    else
    {
        FatalErrorIn("mcParticleCloud::checkMoments()")
            << "Not all moment fields available and no particles present."
            << endl;
    }
}
