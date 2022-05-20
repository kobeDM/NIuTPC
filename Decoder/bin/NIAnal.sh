#!/bin/sh

VERSION=0.0

function usage {
    cat <<EOF
$(basename ${0}) is a tool for ...

Usage:
    $(basename ${0}) [headname] [start_num] [end_num]

Options:
    --version, -v     print $(basename ${0}) version
EOF
}
function version {
    echo "$(basename ${0}) version ${VERSION} "
}
case ${1} in
    help|--help|-h)
        usage
        exit 1
    ;;

    version|--version|-v)
        version
        exit 1
    ;;
esac

if [ $# -ne 3 ]
then
   echo "Wrong input"
   exit
fi
   
# prepare config file
if [ ! -f config.json ]; then
    echo "copy config.json from ${NI_DECODER_DIR}/config/config.json" 
    cp ${NI_DECODER_DIR}/config/config.json ./
fi

#########################################
#  Main
#########################################
COMMAND="NIAnal"

FILE_HEAD=${1}
START_NUM=${2}
END_NUM=${3}

if [ $START_NUM -gt $END_NUM ]
then
    exit
fi

for CURRENT_NUM in `seq -f %03g $START_NUM $END_NUM`
do
    FILE_NAME=${FILE_HEAD}"_"${CURRENT_NUM}".root"
    echo $FILE_NAME
    $COMMAND $FILE_NAME config.json
done

    
    
    

      
