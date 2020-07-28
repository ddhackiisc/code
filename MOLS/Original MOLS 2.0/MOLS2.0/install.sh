echo 'Choose your operating system
If Ubuntu/mint, type 1
If other, type 2'
#If openSuSE, type 2
#If Fedora, type 3
#If other, type 4'
read os

#-----'LOCAL INSTALLATION OF Openbabel2.4.1---------------------------------------
        tar -xzf openbabel-2.4.1.tar.gz
        cd openbabel-2.4.1/
        mkdir ~/MOLS2.0/open-babel-tools
        cmake ../openbabel-2.4.1 -DCMAKE_INSTALL_PREFIX=~/MOLS2.0/open-babel-tools

#echo 'openbabel-2.3.2 installation begins.......'

        make && make install

        cd ../

#-------FPOCKET2 INSTALLATION------------------------------------------------------
	tar -xzf fpocket2.tar.gz

        cd fpocket2

if [ "$os" -eq "1" ]
then
        sed -i 's/$(LINKER) $(LFLAGS) $^ -o $@/$(LINKER) $^ -o $@ $(LFLAGS)/' makefile

fi

        make

        make test

#       make install

        cd ../

#-------chmod-------------------------------------
        chmod 777 lmols
        chmod 777 imolsdock
        chmod 777 reduce
        chmod 777 cpyres
        chmod 777 receptor/centroid
        chmod 777 smi2sdf
        chmod 777 mengine
        chmod 777 interene
#-------------------------------------------------
