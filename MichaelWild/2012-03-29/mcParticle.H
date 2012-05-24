/*---------------------------------------------------------------------------*\
Class
    Foam::mcParticle

Description
    Incomplete example illustrating how to implement a custom particle class.

Author
    Michael Wild

SourceFiles
    mcParticleIO.C
    mcParticle.C

\*---------------------------------------------------------------------------*/

#ifndef mcParticle_H
#define mcParticle_H

#include "particle.H"
#include "contiguous.H"

class mcParticleCloud;

class mcParticle : public particle
{
  // Private Data

    // ORDER IN DESCENDING DATA SIZE!
    // IF YOU UPDATE THIS LIST, ALSO UPDATE mcParticleIO.C!
    scalar m_;   // 8 byte
    vector U_;   // contiguous 8 byte
    vector Ut_;
    scalar rho_;
    /*
    label  n_;   // 4 byte
    bool   b_;   // 1 byte
    */

public:

  TypeName("mcParticle")
  friend class Cloud<mcParticle>;

  //- Persistent storage and helper functions useful during tracking,
  //  e.g. interpolators.
  class trackingData :
    public particle::TrackingData<mcParticleCloud>
  {
      // ...
  };

  // Constructors

    //- Construct from components
    mcParticle(const mcParticleCloud&, const vector& pos, label celli,
               scalar m, const vector& U, scalar rho);
    //- Construct from Istream
    mcParticle(const Cloud<mcParticle>&, Istream&, bool readFields=true);
    //- Clone the particle
    autoPtr<particle> clone() const;
    //- Factory class for parallel transfer
    class iNew
    {
      /* Copy-paste from Foam::particle, s/particle/mcParticle/g */
    };

  // Access

      //- Return the particle mass (statistical weight)
      scalar m() const { return m_; }
      //- Return the particle mass (statistical weight)
      inline scalar& m() { return m_; }
      //- Return the particle velocity
      inline const vector& U() const { return U_; }
      //- Return the particle velocity
      inline vector& U() { return U_; }
      // ...

  // Public Member Functions

    //- Evolve particle for trackTime
    bool move(trackingData&, const scalar trackTime);

  // Patch interactions

    // SKIPPED hit*Patch FUNCTIONS
    hitPatch(/* ...*/); // called before other hit*Patch(/*...*/)

    //- Rotate particle properties (e.g. reflection)
    void transformProperties(const tensor&);
    //- Shift particle properties (e.g. cyclic patch)
    void transformProperties(const vector&);

  // I/O

    //- Read properties from files <time>/lagrangian/<cloudName>/<prop>
    static void readFields(Cloud<mcParticle>&);
    //- Write properties to files <time>/lagrangian/<cloudName>/<prop>
    static void writeFields(const Cloud<mcParticle>&);
    //- Write properties to stream (e.g. for parallel transfer)
    friend Ostream& operator<<(Ostream&, const mcParticle&);
};

// IF CONTAINS NON-CONTIGUOUS DATA
/*
template<>
inline bool contiguous<mcParticle>()
{
  return false;
}
*/

#endif
