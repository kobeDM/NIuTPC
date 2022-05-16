DIR=`pwd`
echo "NIuTPC analysis package in ${DIR} is set"
export NI_DECODER_DIR=${DIR}
export PATH=${DIR}/bin:${PATH}
export PATH=${DIR}/macros:${PATH}
