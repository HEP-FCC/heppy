

source /cvmfs/fcc.cern.ch/sw/latest/setup.sh

export HEPPY=$PWD/heppy
export PATH=$HEPPY/bin:$PATH
export PYTHONPATH=$PWD:$PYTHONPATH

# need this for heppy's context discovery. TODO: get rid of context discovery in heppy
export FCCEDM="unused"
export PODIO="unused"
export FCCPHYSICS="unused"

