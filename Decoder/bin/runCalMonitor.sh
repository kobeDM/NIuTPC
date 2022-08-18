#!/bin/bash

INPUTLIST=$1
OUTPUTDIR=$2

# prepare rootlogon
if [ ! -f rootlogon.C ]; then
    echo "copy rootlogon.C from ${NI_DECODER_DIR}/macros/rootlogon.C" 
    cp ${NI_DECODER_DIR}/macros/rootlogon.C ./
fi

# prepare config file
if [ ! -f config.json ]; then
    echo "copy config.json from ${NI_DECODER_DIR}/config/config.json" 
    cp ${NI_DECODER_DIR}/config/config.json ./
fi

ls -lhv

# execute root macro
COMMAND=${NI_DECODER_DIR}/macros/calMonitor.cc
root -q -b -l "${COMMAND}(\"${INPUTLIST}\",\"config.json\",\"${OUTPUTDIR}\")"

