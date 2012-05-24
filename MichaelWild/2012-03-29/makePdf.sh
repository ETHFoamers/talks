#!/bin/sh

set -e

enscript -Ma4 -U2 -Ecpp -hj --margins=15:15:15:15 -o mcParticle.ps \
	mcParticle.H mcParticle.C mcParticleIO.C \
	mcParticleCloud.H mcParticleCloud.C

ps2pdf12 mcParticle.ps

rm -f mcParticle.ps
