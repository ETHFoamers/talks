#include "mcParticle.H"
#include "IOstreams.H"
#include "mcParticleCloud.H"

Foam::mcParticle::mcParticle
(
    const Cloud<mcParticle>& cloud,
    Istream& is,
    bool readFields
)
:
    particle(cloud.pMesh(), is, readFields)
{
    if (readFields)
    {
        if (is.format() == IOstream::ASCII)
        {
            m_ = readScalar(is);
            is  >> U_
                >> Ut_
                >> rho_
                ;
        }
        else
        {
            is.read
            (
                reinterpret_cast<char*>(&m_),
                sizeof(m_) + sizeof(U_) + sizeof(Ut_) + sizeof(rho_)
            );
        }
    }
    // Check state of Istream
    is.check("mcParticle::mcParticle(Istream&)");
}


void Foam::mcParticle::readFields(Cloud<mcParticle>& c)
{
    if (!c.size())
    {
        return;
    }
    particle::readFields(c);

    mcParticleCloud& mcpc = refCast<mcParticleCloud>(c);

    IOField<scalar> m(c.fieldIOobject("m", IOobject::MUST_READ));
    c.checkFieldIOobject(c, m);

    IOField<vector> U(c.fieldIOobject("U", IOobject::MUST_READ));
    c.checkFieldIOobject(c, U);

    IOField<scalar> rho(c.fieldIOobject("rho", IOobject::MUST_READ));
    c.checkFieldIOobject(c, rho);

    label i = 0;
    forAllIter(Cloud<mcParticle>, c, iter)
    {
        mcParticle& p = iter();

        p.m_ = m[i];
        p.U_ = U[i];
        p.rho_ = rho[i];
        ++i;
    }
}


void Foam::mcParticle::writeFields(const Cloud<mcParticle>& c)
{
    particle::writeFields(c);

    const mcParticleCloud& mcpc = refCast<const mcParticleCloud>(c);

    label np = c.size();

    IOField<scalar> m(c.fieldIOobject("m", IOobject::NO_READ), np);
    IOField<vector> U(c.fieldIOobject("U", IOobject::NO_READ), np);
    IOField<scalar> rho(c.fieldIOobject("rho", IOobject::NO_READ), np);

    label i = 0;
    forAllConstIter(Cloud<mcParticle>, c, iter)
    {
        const mcParticle& p = iter();

        m[i] = p.m_;
        U[i] = p.U_;
        rho[i] = p.rho_;
        i++;
    }

    m.write();
    U.write();
    rho.write();
}


Foam::Ostream& Foam::operator<<(Ostream& os, const mcParticle& p)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const particle&>(p)
            << token::SPACE << p.m_
            << token::SPACE << p.U_
            << token::SPACE << p.Ut_
            << token::SPACE << p.rho_;
    }
    else
    {
        os  << static_cast<const particle&>(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.m_),
            sizeof(p.m_) + sizeof(p.U_) + sizeof(p.Ut_) + sizeof(p.rho_)
        );
    }
    // Check state of Ostream
    os.check("Ostream& operator<<(Ostream&, const mcParticle&)");
    return os;
}
