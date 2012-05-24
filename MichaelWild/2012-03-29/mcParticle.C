#include "mcParticle.H"
#include "mcParticleCloud.H"

namespace Foam
{
    defineTypeNameAndDebug(mcParticle, 0);
}


Foam::mcParticle::mcParticle
(const mcParticleCloud& c, const vector& pos,
label celli, scalar m, const vector& U, scalar rho);
:
    particle(c.pMesh(), pos, celli),
    m_(m),
    U_(U),
    Ut_(),
    rho_(rho),
    n_(0),
    b_(false)
{}


bool Foam::mcParticle::move
(
    mcParticle::trackingData& td,
    const scalar trackTime
)
{
    td.switchProcessor = false;
    td.keepParticle = true;
    const mcParticleCloud& mcpc = refCast<mcParticleCloud>(td.cloud());
    const polyMesh& mesh = mcpc.pMesh();
    const polyBoundaryMesh& pbMesh = mesh.boundaryMesh();

    scalar tEnd = (1.0 - stepFraction())*trackTime;
    scalar dtMax = tEnd;

    // At beginning of time step, update velocity
    if (stepFraction() < SMALL)
    {
        /* SKIPPED U_ +=  .... ; */
    }

    // Compute tracking velocity (i.e. handle 2D and wedge cases)
    Ut_ = U_;
    meshTools::constrainDirection(mesh, mesh.solutionD(), Ut_);
    point destPos = position() + tEnd * Ut_;
    if (mcpc.isAxiSymmetric())
    {
        vector rotatedCentreNormal = mcpc.axis()^destPos;
        rotatedCentreNormal /= mag(rotatedCentreNormal);
        tensor T = rotationTensor(rotatedCentreNormal, mcpc.centrePlaneNormal());
        transformProperties(T);
        destPos = transform(T, destPos);
        // constrain to kill numerical artifacts
        meshTools::constrainDirection(mesh, mesh.geometricD(), destPos);
        Ut_ = (destPos - position())/tEnd;
    }

    while (td.keepParticle && !td.switchProcessor && tEnd > SMALL)
    {
        scalar dt = min(dtMax, tEnd);
        destPos = position() + dt*Ut_;
        // do actual tracking
        scalar tf = trackToFace(destPos, td);
        ++nSteps_;
        dt *= tf;
        tEnd -= dt;
        stepFraction() = 1.0 - tEnd/trackTime;

        if (onBoundary() && td.keepParticle)
        {
            if (isA<processorPolyPatch>(pbMesh[patch(face())]))
            {
                td.switchProcessor = true;
            }
        }
    }
    return td.keepParticle;
}


void Foam::mcParticle::transformProperties (const tensor& T)
{
    particle::transformProperties(T);
    U_ = transform(T, U_);
    Ut_ = transform(T, Ut_);
}


void Foam::mcParticle::transformProperties(const vector& separation)
{
    particle::transformProperties(separation);
}

// SKIPPED hit*Patch(...) FUNCTIONS
