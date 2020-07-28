#written by Sam Paul D. 
#script to find the centroid of the first five pockets of fpocket.
#this script is used by 'auto' option of  gridcentre
#date : 20th November 2012

#rm receptor.pdb
./cpyres #copy receptor(protein) to 'receptor.pdb'
mv receptor.pdb receptor/

cd receptor

rm -rf *_out

fpocket -f receptor.pdb
echo 'fpocket done'

#copy "ATOM" section from pocket no.1 to 'pocket.pdb' 
#'centroid' is the executable of prelims.f, which is set to always find the centroid of 'pocket.pdb'
 
grep "ATOM" *_out/pockets/pocket0_atm.pdb > pocket.pdb
./centroid > pockcentroid1
rm part.pdb
#rm pocket.pdb
#grep "ATOM" *_out/pockets/pocket1_atm.pdb > pocket.pdb
#./centroid > pockcentroid2
#rm part.pdb
#rm pocket.pdb
#grep "ATOM" *_out/pockets/pocket2_atm.pdb > pocket.pdb
#./centroid > pockcentroid3
#rm part.pdb
#rm pocket.pdb
#grep "ATOM" *_out/pockets/pocket3_atm.pdb > pocket.pdb
#./centroid > pockcentroid4
#rm part.pdb
#rm pocket.pdb
#grep "ATOM" *_out/pockets/pocket4_atm.pdb > pocket.pdb
#./centroid > pockcentroid5
#rm part.pdb
#rm pocket.pdb

echo 'pocket copied and centroid calculated'

#values of x,y,z of the centroid of the first five pockets are copied to file 'centroid_pocket_*'
sh textbox_pocket1_1.sh > centroid_pocket_1_1
sh textbox_pocket1_2.sh > centroid_pocket_1_2
sh textbox_pocket1_3.sh > centroid_pocket_1_3

#sh textbox_pocket2_1.sh > centroid_pocket_2_1
#sh textbox_pocket2_2.sh > centroid_pocket_2_2
#sh textbox_pocket2_3.sh > centroid_pocket_2_3

#sh textbox_pocket3_1.sh > centroid_pocket_3_1
#sh textbox_pocket3_2.sh > centroid_pocket_3_2
#sh textbox_pocket3_3.sh > centroid_pocket_3_3

#sh textbox_pocket4_1.sh > centroid_pocket_4_1
#sh textbox_pocket4_2.sh > centroid_pocket_4_2
#sh textbox_pocket4_3.sh > centroid_pocket_4_3

#sh textbox_pocket5_1.sh > centroid_pocket_5_1
#sh textbox_pocket5_2.sh > centroid_pocket_5_2
#sh textbox_pocket5_3.sh > centroid_pocket_5_3

echo 'centroid copied'

#rm pockcentroid1
#rm pockcentroid2
#rm pockcentroid3
#rm pockcentroid4
#rm pockcentroid5



