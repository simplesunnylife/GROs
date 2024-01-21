# GROs
code and documents for GRO projects

# Description：

Demo code for VLDB paper :

Title

The code is developed based library GUNDAM https://github.com/MinovskySociety/GUNDAM

Main repo for graph rule discovery is being refactored and optimized, which will be released soon. https://github.com/MinovskySociety/GraphRules

Install dependencies on Ubuntu
GCC version: 7.4.0 or above, support of c++17 standard required.
```bash
Install mpi:

sudo apt-get install openmpi-bin openmpi-doc libopenmpi-dev
Install glog:

sudo apt-get install libgoogle-glog-dev
Install gflags:

sudo apt-get install libgflags-dev
Install yaml:

sudo apt-get install libyaml-cpp-dev
Compile
mkdir build && cd ./build
cmake ../
make all -j
mv ./dataset/* ./build.

```
