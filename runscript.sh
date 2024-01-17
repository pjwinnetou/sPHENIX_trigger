#! /bin/bash
export USER="$(id -u -n)"
export LOGNAME=${USER}
export HOME=/sphenix/user/${LOGNAME}

source /opt/sphenix/core/bin/sphenix_setup.sh -n new
SEG=${1}
echo $SEG
./calcEff $SEG

echo all done
