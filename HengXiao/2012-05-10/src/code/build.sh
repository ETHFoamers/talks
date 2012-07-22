#!/bin/sh

set -e

enscript -U2 -Ecpp -hj --margins=15:15:15:15 -o durbin-code.ps \
    Durbin.C muGenElliptic.C

ps2pdf12 durbin-code.ps

rm -f durbin-code.ps
