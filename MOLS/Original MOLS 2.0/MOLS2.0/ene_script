obenergy -ff MMFF94 cal_energy.mol2 > ene
grep 'TORSIONAL\|VAN DER WAALS\|ELECTROSTATIC' ene | tail -n 3 >babelene.txt
awk '/TORSIONAL/ {print $5}' babelene.txt  > output1
awk '/VAN/ {print $7}' babelene.txt >> output1
awk '/ELECTROSTATIC/ {print $5}' babelene.txt >> output1


