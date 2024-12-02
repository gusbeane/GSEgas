#!/bin/bash

ICS_BASE=/n/holylfs05/LABS/hernquist_lab/Users/abeane/GSEgas/ics/MW6iso/lvl5-dens-grid
COMMIT=763e6294d
BRANCH=SMUGGLEdev-dt_AB
CONFIG=Config.sh
LVL=5

if (( $# != 3 )); then
    echo -e "Need beta, RC, vphi" >&2
    exit 1  # Return with an error status
fi

BETA=$1
RC=$2
VPHI=$3

source load-modules.sh

if test -d lvl${LVL}-beta${BETA}-RC${RC}-vphi${VPHI}; then
  echo "Directory already exists."
  exit 1
fi

mkdir lvl${LVL}-beta${BETA}-RC${RC}-vphi${VPHI}
cd lvl${LVL}-beta${BETA}-RC${RC}-vphi${VPHI}

    cp ../load-modules.sh .
    ln -s ${ICS_BASE}/beta${BETA}/RC${RC}-vphi${VPHI} ICs

    cp ../param/param_lvl${LVL}.txt ./
    sed -i "s/MW/MW-beta${BETA}-RC${RC}-vphi${VPHI}/g" param_lvl${LVL}.txt
    
    cp ../$CONFIG Config.sh

    git clone --depth=32 --branch=${BRANCH} git@bitbucket.org:volkerspringel/arepo.git

    cd arepo
        git checkout $COMMIT || exit 1
        ln -s ../Config.sh ./
    cd ../

    for jn in casc; do
        cp ../job_lvlN-${jn}.sh job_lvl${LVL}-${jn}.sh
        sed -i "s/MW6btBTrcRCvpVP/MW6bt${BETA}rc${RC}vp${VPHI}/g" job_lvl5-${jn}.sh
        sed -i "s/lN/l${LVL}/g" job_lvl${LVL}-${jn}.sh
    done

cd ../

