#fname=inp.mol
#for x in `cat $fname`
#do
#babel -imol $x -osmi inp.smi
#done

babel -imol ligand.mol -osmi inp.smi
