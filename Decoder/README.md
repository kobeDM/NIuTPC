---
# Decoder0.2_cmake
## author	T.Shimada
## date	2020.04.02
---

cmake v3.3

##How to compile

If there is not "build" directry, you should make directory.

```
$ pwd
/path/to/Decoder0.2_cmake

$ ls
Makefile README.md cmake/ source/

$ make clean
$ make Event

$ mkdir bin
$ mkdir build
$ cd build
$ cmake ../source
$ make
$ make install
$ ls
Makefile README.md bin/ build/ cmake/ source/

$ cp source/Event/EventDict_rdict.pcm bin/
```
After compiling, write .bashrc
```
export PATH=/path/to/Decoder/bin:$PATH


