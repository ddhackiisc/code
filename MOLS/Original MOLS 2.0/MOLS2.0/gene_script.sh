
obenergy -ff GAFF gcal_energy.mol2 > gene

#grep 'BOND STRETCHING ENERGY' ene > babelene.txt
#grep 'ANGLE BENDING ENERGY' ene >> babelene.txt
#grep 'TORSIONAL ENERGY' ene >> babelene.txt
#grep 'IMPROPER-TORSIONAL ENERGY' ene >> babelene.txt
#grep 'VAN DER WAALS ENERGY' ene >> babelene.txt
#grep 'ELECTROSTATIC ENERGY' ene >> babelene.txt
#grep 'TOTAL ENERGY' ene >> babelene.txt

grep 'TORSIONAL\|VAN DER WAALS\|ELECTROSTATIC' gene | tail -n 4 >babelene.txt
awk '/TORSIONAL/ {print $5}' babelene.txt  > output1
awk '/VAN/ {print $7}' babelene.txt >> output1
awk '/ELECTROSTATIC/ {print $5}' babelene.txt >> output1


#awk '/TOTAL ENERGY/ {print $4}' gene > output1


