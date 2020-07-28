#echo 'type the location of the MOLS2.0 executable
#e.g. /home/user-name/MOLS2.0'
#read dir
dir=$(pwd)

echo "[Desktop Entry]
Version=2.0
Type=Application
Terminal=false
Exec=java -jar $dir/mols.jar
Name=MOLS2.0
Path=$dir
Comment=a protein-ligand docking algorithm
#Icon=/usr/share/new.png            #you can add an icon
Icon= $dir/mols_icon.png
Categories=Education" > mols.desktop
cp mols.desktop ~/.local/share/applications
exit

