/*---------------------------------------------------------------------------*\
Class
    Foam::mcParticleCloud

Description
    Incomplete example illustrating how to implement a custom particle cloud.

Author
    Michael Wild

SourceFiles
    mcParticleCloud.C

\*---------------------------------------------------------------------------*/

#ifndef mcParticleCloud_H
#define mcParticleCloud_H

#include "Cloud.H"
#include "mcParticle.H"
#include "dictionary.H"

class mcParticleCloud:
    public Cloud<mcParticle>
{
  // Private Data

    //- FV velocity
    const volVectorField& Ufv_;
    //- FV pressure
    const volScalarField& pfv_;
    //- Averaged 0th statistical moment
    DimensionedField<scalar, volMesh> mMom_;
    //- Averaged 1st statistical moment of inverse density (volume)
    DimensionedField<vector, volMesh> VMom_;
    //- Averaged 1st statistical moment of velocity
    DimensionedField<vector, volMesh> UMom_;
    //- Averaged 2nd statistical moment of fluctuating velocity
    DimensionedField<symmTensor, volMesh> uuMom_;
    //- Extracted, time-averaged density
    volScalarField rho_;
    //- Extracted, time-averaged particle-mass density
    volScalarField pmd_;
    //- Extracted, time-averaged velocity
    volVectorField U_;
    //- Extracted, time-averaged turbulent stress tensor
    volSymmTensorField Tau_;
    //- Extracted, time-averaged TKE
    volScalarField k_;

  // Private Member Functions

    //- Update moments and the quantities remembered by particles
    void updateCloudPDF(scalar existWt);

    //- Initialize statistical moments
    void initMoments();

    //- Disallow default bitwise copy construct
    mcParticleCloud(const mcParticleCloud&);
    //- Disallow default bitwise assignment
    void operator=(const mcParticleCloud&);

public:

  // Constructors

    //- Construct from components
    // If U, p or rho are NULL, they are looked up by the names defined by the
    // UName, pName rhoName entries (defaulting to U, p and rho, respectively)
    mcParticleCloud
    (
      const fvMesh& mesh,
      const dictionary& dict,
      const word& cloudName = "defaultCloud",
      const volVectorField* U = 0,
      const volScalarField* p = 0,
      const volScalarField* rho = 0
    );

  // Access

    //- Returns the FV TKE field (SKIPPED HERE)
    const volScalarField& kfv() const;

    //- Returns the FV turbulence frequency field (SKIPPED HERE)
    const volScalarField& omega() const;

    //- Returns true if the case uses a wedge geometry
    bool isAxiSymmetric() const;

    //- Axis of the wedge (only defined if isAxiSymmetric())
    const vector& axis() const;

    //- Centre plane of the wedge (only defined if isAxiSymmetric())
    const vector& centrePlaneNormal() const;

  // Public Member Functions

    //- Initialize the particle cloud
    void initReleaseParticles();

    //- Randomly generate N particles in celli, with provided cell-based
    // values and the scale of velocity fluctuation
    void particleGenInCell
    (
      label celli,
      label N,
      scalar m,
      const vector& U,
      const vector& urms
    );

    //- Evolve the particles
    void evolve();
};

#endif
