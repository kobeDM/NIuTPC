DIR=`pwd`
echo ${DIR}

ls -lh

echo "initial setup..."
make Event

mkdir build bin
cd build
cmake ../source
make && make install

ls -lh

cd ${DIR} 
cp source/Event/EventDict_rdict.pcm bin/

export NI_DECODER_DIR=${DIR}
export PATH=${DIR}/bin:${PATH}
