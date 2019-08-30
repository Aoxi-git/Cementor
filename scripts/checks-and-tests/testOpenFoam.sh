#! /bin/bash

ls -la /root/OpenFOAM/OpenFOAM-6/etc/bashrc
source  /root/OpenFOAM/OpenFOAM-6/etc/bashrc

cd /builds/yade-dev/
rm -rf Yade-OpenFOAM-coupling
git clone https://github.com/dpkn31/Yade-OpenFOAM-coupling.git
cd Yade-OpenFOAM-coupling
git checkout yadetest
./Allclean
./Allwmake

cd /builds/yade-dev/trunk/examples/openfoam/example_icoFoamYade
ln -s /builds/yade-dev/trunk/install/bin/yade-ci ./yadeimport.py
blockMesh
decomposePar
mkdir yadep

mpiexec --allow-run-as-root -n 1 python3 scriptYade.py : -n 2 icoFoamYade -parallel

