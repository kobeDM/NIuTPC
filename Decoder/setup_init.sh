DIR=`pwd`
echo ${DIR}

ls -lh

echo "initial setup..."
make Event

mkdir build
cd build
cmake ../source
make && make install

ls -lh

cd ${DIR} 
cp source/Event/EventDict_rdict.pcm bin/

export PATH=${DIR}/bin:${PATH}
