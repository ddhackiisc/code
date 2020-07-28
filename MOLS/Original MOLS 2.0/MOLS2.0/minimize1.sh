obminimize -n 2500 -ff MMFF94 -cg mols.mol2 >min.pdb
#babel -ipdb min.pdb -omol2 min.mol2

#obenergy -ff MMFF94 min.mol2 > ene11
#grep 'TORSIONAL\|VAN DER WAALS\|ELECTROSTATIC' ene11 | tail -n 3 >babelene.txt
#awk '/TORSIONAL/ {print $5}' babelene.txt  > output2
#awk '/VAN/ {print $7}' babelene.txt >> output2
#awk '/ELECTROSTATIC/ {print $5}' babelene.txt >> output2

