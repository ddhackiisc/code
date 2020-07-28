#echo 'type the location of the Jmolsdock executable
#e.g. /home/programs/Jmolsdock'
#read dir
dir=$(pwd)
echo 'Choose your operating system
If Ubuntu/mint, type 1
If openSuSE, type 2
If Fedora, type 3
If other, type 4'
read os

#-------zlib installation-----------------------------
if [ "$os" -eq "1" ] || [ "$os" -eq "2" ] || [ "$os" -eq "3" ]
then
	tar -xzf zlib.tar.gz

	cd zlib
	
	./configure

	make test

	sudo make install

	cd ../
fi


#echo '**********************************************'
#echo 'GLOBAL Installation of Openbabel2.2.0
#echo '**********************************************'
tar -xzf openbabel-2.2.0.tar.gz

	cd openbabel-2.2.0/

#echo 'openbabel-2.2.0 installation begins.......'

	./configure 

	sudo make

#echo 'sudo install openbabel'
#echo 'enter you password'	
	sudo make install
	cd ../


if [ "$os" -eq "1" ] || [ "$os" -eq "2" ]
then
        sudo /sbin/ldconfig
fi




#	echo '******************************************'
#       echo 'openbabel-2.2.0 installation successful!!!'
#	echo '******************************************'
	tar -xzf fpocket2.tar.gz

	cd fpocket2

if [ "$os" -eq "1" ]
then 
	sed -i 's/$(LINKER) $(LFLAGS) $^ -o $@/$(LINKER) $^ -o $@ $(LFLAGS)/' makefile

fi

        sudo make

	make test

	sudo make install

	cd ../
	#echo '**********************************'
        #echo 'fpocket installation successful!!!'
	#echo '**********************************'
#else
echo "[Desktop Entry]
Version=2.0
Type=Application
Terminal=true
Exec=java -jar $dir/Jmolsdock.jar
Name=Jmolsdock
Path=$dir
Comment=a protein-ligand docking algorithm
#Icon=/usr/share/new.png            #you can add an icon
Categories=Education" > Jmolsdock.desktop
sudo cp Jmolsdock.desktop /usr/share/applications
#	exit

