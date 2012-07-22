Class
    Foam::incompressible::RASModels::muGenElliptic

Description
    An abstract virtual class for elliptic relaxation RSTM models.
    Examples include:

    Durbin Elliptic Relaxation Model.
    
    Durbin, A Reynolds stress model for near-wall turbulence, J. Fluid
    Mech. vol 249, 1993.

    and 
    Hanjalic Elliptic Blending Model

    Manceau and Hanjalic, Elliptic blending model: a new
    near-wall Reynolds-stress turbulence closure. Physics of Fluids 14
    (74), 2002.
    
    The default model coefficients correspond to the following:
SourceFiles
    muGenElliptic.C

\*---------------------------------------------------------------------------*/

#ifndef muGenElliptic_H
#define muGenElliptic_H

#include "muRASModel.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class muGenElliptic Declaration
\*---------------------------------------------------------------------------*/

class muGenElliptic
:
    virtual public muRASModel
{
    // Private member functions:

    // Disallow default bitwise copy construct and assignment
    muGenElliptic(const muGenElliptic&); 
    muGenElliptic& operator=(const muGenElliptic&); 

    //Protected data

    protected:
        // Model coefficients
            const fvMesh& mesh_;
            dimensionedScalar CmuKE_;

            // LRR-IP model related coeffs
            dimensionedScalar Clrr1_;
            dimensionedScalar Clrr2_;

            dimensionedScalar C1_;
            dimensionedScalar C2_;
            dimensionedScalar Cmu_;
            dimensionedScalar sigmaEps_;
            dimensionedScalar sigmaK_;

            // coeffs related to elliptic relaxation
            dimensionedScalar CL_;
            dimensionedScalar CEta_;
            dimensionedScalar yStarLim_; 
            const Switch      implicitDiv_;

            volScalarField KolmogorovFlag_;
    
      // Fields
            volSymmTensorField R_;
            volScalarField k_;
            volScalarField epsilon_;

       // Forcing fields
            volSymmTensorField QR_;
            volScalarField     Qepsilon_;
            Switch          directImpose_;
            labelList       lesCellLables_;
            scalarField     lesCellValuesK_;
            scalarField     lesCellValuesEpsilon_;
            symmTensorField lesCellValuesR_;
            label           imposeTurbEvery_;

        // Private member functions
        tmp<volTensorField> symm2full(volSymmTensorField& symm ) const;
        void updateKolmogorovFlag();
        tmp<volScalarField> T() const;
        tmp<volScalarField> L() const;

public:

    // Constructors

        //- Construct from components
        muGenElliptic
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport,
            const volSymmTensorField & RAvg,
            const volScalarField & epsilonAvg
        );


    //- Destructor
    virtual ~muGenElliptic()
    {}


    // Member Functions

        //- Return the turbulence viscosity
        virtual tmp<volScalarField> nut() const
        {
            return CmuKE_*sqr(k_)/(epsilon_ + epsilonSmall_);;
        }

        //- Return the effective stress tensor including the laminar stress
        virtual tmp<volSymmTensorField> devReff() const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevReff(volVectorField& U) const;

        //- update Kolmogorov flag field only
        virtual void correct();

        //- Read RASProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

