#!/bin/sh

set -e

enscript -E --color -whtml --toc -p mcParticle.html \
	mcParticle.H mcParticle.C mcParticleIO.C \
	mcParticleCloud.H mcParticleCloud.C
