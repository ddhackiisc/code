#-------cmake installation----------------------------
        mkdir ~/MOLS2.0/cmake_mols
        tar -zxvf cmake-3.5.1.tar.gz
        cd cmake-3.5.1/
        ./bootstrap
        make
        cd ../
#-------zlib installation-----------------------------
        mkdir ~/MOLS2.0/zlib_mols
        tar -xzf zlib.tar.gz
        cd zlib
        ./configure --prefix=~/MOLS2.0/zlib_mols
        make test
        make install
        cd ../
#-----------------------------------------------------
