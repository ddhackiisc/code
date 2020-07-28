#!/usr/bin/perl -w
#
# Script which determines the rotatable bonds in mol2 file(s),
# written by Visvaldas Kairys <kairys (at) uma.pt>
# Created in about 2000-2001, when I was a postdoc at 
# Center for Advanced Research in Biotechnology / University of Maryland
#
# The criteria for prohibited flexible bonds are taken from
# S. Makino, I. D. Kuntz  J. Comput. Chem., vol. 18, p. 1812-1825 (1997).
#
# In other words, the script pays attention to the atom types: does not rotate
# double bonds, etc.
# Also finds terminal bonds and considers them non-rotatable (example: C-H bond).
# In addition, considers methyls as non-rotatable.
# After little editing, some of the restrictions could be removed.
# If you want this behavior to be modified, drop me a line, 
# The scheme is not perfect (for example, -NH3 groups are rotatable,
# actually this could be added to the script, but I didn't get to that), 
# but hopefully is usable.
# 
# Of course, bonds inside the rings are not rotatable.
#
# Ring detection algorithm taken from:
# Th. Hanser, Ph. Jauffret, G. Kaufmann J. Chem. Inf. Comput. Sci. 1996, vol 36, p. 1146-1152.

# Appreciation to Mike Potter for tips and Mike Gilson for guidance.

#
# Usage: findrotatable.pl file1.mol2 <file2.mol2 file3.mol2 ...>
#

# The info about rotatable bonds is output to the screen, and also to 
# file(s) with extension *.aux, one for each molecule.


# Copyright (C)  2006  Visvaldas Kairys

#This program is free software; you can redistribute it and/or
#modify it under the terms of the GNU General Public License
#as published by the Free Software Foundation; either version 2
#of the License, or (at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program; if not, write to the Free Software
#Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

if (@ARGV < 1 ) {
    die "Usage: $0 mol2file(s)\n";
}

foreach $file (@ARGV){
	next unless ($file =~ /mol2/i);
        $file =~ s/\..*?$//; # remove file extension
        open(MOL2FILE,"$file.mol2") or die "Cannot open file $file.mol2 $!\n";
        print ">>>>>>Processing file $file.mol2...\n";
        my %types=();
        my %spectypes=();
        my @m_bonds=();
        my @ring=();
        my @ringbonds=();
	$tmplineno=0;
	$firstatomline=0;
	$lastatomline=0;
	$firstbondline=0;
	$lastbondline=0;
        while(<MOL2FILE>){
		if(/TRIPOS>MOLEC/){
			$tmplineno=$.;
		}
		if($.==$tmplineno+2){
                       	chomp;
                       	@tmp=split;
			$natoms=$tmp[0];
			$nbonds=$tmp[1];
			print "Number of atoms: $natoms, bonds: $nbonds\n";
		}
		if(/TRIPOS>ATOM/){
                        $firstatomline=$.+1;
                        $lastatomline=$.+$natoms;
                        #print "first $firstatomline last $lastatomline\n";
                        next;
                }
                if($. == $firstatomline .. $. == $lastatomline){
                        #print "linenumber $., \$lastatomline= $lastatomline, \$natoms=$natoms\n";
                        chomp;
                        @tmp=split;
                        #create a list containing atom number as a key and type as a value
                        $types{$tmp[0]}=$tmp[5];
                }
		if(/TRIPOS>BOND/){
			$firstbondline=$.+1;
			$lastbondline=$.+$nbonds;
			#print "first $firstbondline last $lastbondline\n";
			next;
		}

		if($. == $firstbondline .. $. == $lastbondline){
			#print "linenumber $., \$lastbondline= $lastbondline, \$nbonds=$nbonds\n";
                        chomp;
                        @tmp=split;
                        #create an array of arrays containing bond data
                        push @m_bonds, [ @tmp ];
			#print "mbond $tmp[0] $tmp[1] $tmp[2]\n";
		}
	}

	#foreach $i (keys %types){
	#	print "types $i $types{$i}\n";
        #        }

	#print "items in \@m_bonds= $#m_bonds\n";
	#foreach $i (0 .. $#m_bonds){
	#	print "mbonds $i $m_bonds[$i]->[0] $m_bonds[$i]->[1] $m_bonds[$i]->[2]\n";
        #        }

        
        	my %neighb=();
			foreach $ii (keys %neighb){
				print " init $ii : @{$neighb{$ii}} \n";
			}
	foreach $i (keys %types){
		for $j (0 .. $#m_bonds){  #find all bonds connected to that carbon
			$i1 = $m_bonds[$j]->[1];
			$i2 = $m_bonds[$j]->[2];
			if($i1 == $i || $i2 == $i){
				#print "$i1 $i2 $m_bonds[$j]->[3]\n";
                                $atom=$i2 if($i1 == $i);
                                $atom=$i1 if($i2 == $i);
				push @{$neighb{$i}},$atom;
				}
		}
	}
# before we start looking for non-rotatable bonds we determine if some atoms belong to some more complicated atom types
LSPC:	foreach $i (keys %types){
		$count = 0;
		my @values=@{$neighb{$i}};
			for $j (0 .. $#values){
#uncomment for debug#				print " $i neigh types: $types{$values[$j]}\n";
			}
#               planar amine
		if($types{$i} =~ /N\.pl3|N\.2/){
			# there have to be exactly 2 hydrogens connected to it 
			for $j (0 .. $#values){
				$count += 1 if $types{$values[$j]} =~ /H/;
			}
			if($count==2){
				$spectypes{$i}='planar_nh2';
				next LSPC;
			}
		}

#               methyl
#               and as a special case methyl connected to N,O or S.
		if($types{$i} =~ /C\.3/){
			# there have to be exactly 3 hydrogens connected to it 
			my $noscount=0;
			for $j (0 .. $#values){
				$count += 1 if $types{$values[$j]} =~ /H/;
				$noscount += 1 if $types{$values[$j]} =~ /N\.|O\.|S\./;
			}
			if($count==3){
				if ($noscount==1){
					$spectypes{$i}='nos_methyl';
				}else{
					$spectypes{$i}='methyl';
				}
				next LSPC;
			}
		}

#               amidine
		if($types{$i} =~ /C\.2|C\.cat/){
			for $j (0 .. $#values){
				$count += 1 if $types{$values[$j]} =~ /N.pl3|N.2/;
			}
			if($count==2){
				$spectypes{$i}='amidine';
				next LSPC;
			}
		}

#               sulfonic acid
		if($types{$i} =~ /S\.3/){
			for $j (0 .. $#values){
				$count += 1 if $types{$values[$j]} =~ /O\.co2/;
			}
			if($count==3){
				$spectypes{$i}='sulfonate';
				next LSPC;
			}
		}

#               nitro,nitroso
		if($types{$i} =~ /N\..*/){
			for $j (0 .. $#values){
				$count += 1 if $types{$values[$j]} =~ /O\.2/;
			}
			if($count>=1){
				$spectypes{$i}='nitro';
				next LSPC;
			}
		}

#               carboxylic acid
		if($types{$i} =~ /C\.2|C\.cat/){
			for $j (0 .. $#values){
				$count += 1 if $types{$values[$j]} =~ /O\.co2/;
			}
			if($count==2){
				$spectypes{$i}='carboxy';
				next LSPC;
			}
		}

		$spectypes{$i}='none';
	}
#uncomment for debug#	foreach $i (keys %spectypes){
#uncomment for debug#		print "spectype for $i $spectypes{$i}\n" unless ($spectypes{$i} =~ 'none');
#uncomment for debug#                }

#       initialize edges
	my %edge=();
	for $m (0 .. $#m_bonds){ 
		$edge{$m}->[0] = $m_bonds[$m]->[1];
		$edge{$m}->[1] = $m_bonds[$m]->[2];
#uncomment for debug#			print "Edges $m ------ @{$edge{$m}} \n";
	}

#       initialize vertices; use the number of connectivities (neighbors) as the value.
	my %vertex=();
	foreach $m (keys %types){ 
		$nb=@{$neighb{$m}} ;
		$vertex{$m}=$nb;
#uncomment for debug#		print "Vertices $m ======= $vertex{$m} \n";
	}

	$nv=scalar (keys %vertex);
	$ne=scalar (keys %edge);
#uncomment for debug#	print "Renewed vertices: $nv edges : $ne \n";

	$oneneighbor=1;
	while($oneneighbor>0){
#       delete appendages to the rings. Update %vertex and %edge hashes
		foreach $m (keys %vertex){ 
			if($vertex{$m}==1){
#uncomment for debug#				print "entry $m deleted in vertex\n";
				delete $vertex{$m};
				foreach  $n (keys %edge){ 
					$k1=$edge{$n}->[0];
					$k2=$edge{$n}->[1];
					if($k1==$m || $k2==$m){
						delete $edge{$n};
#uncomment for debug#						print "entry $n deleted in edge\n";
						if($k1==$m){
							$vertex{$k2} -= 1
						}
						if($k2==$m){
							$vertex{$k1} -= 1
						}
					}
				}
			}
		}
		$oneneighbor=0;
		foreach $m (keys %vertex){ 
#uncomment for debug#			print "new vertex: $m $vertex{$m}\n";
			$oneneighbor += 1 if $vertex{$m}==1;
		}
		$nv=scalar (keys %vertex);
		$ne=scalar (keys %edge);
#uncomment for debug#		print "Renewed vertices: $nv, edges : $ne, atoms with 1 neighbor : $oneneighbor\n";
	}
	while(scalar (keys %vertex) > 0){
	my @label=();
#       reduce the graph further by collapsing other vertices

#       sort the vertex hash by increasing number of neighbors; pick just one, first vertex
	foreach $m (sort {$vertex{$a} <=> $vertex{$b} } keys %vertex) {
		#print "sorted: $m ---> $vertex{$m}\n";
		$curvert=$m;
		last;
	}
#uncomment for debug#	print "vertex $curvert picked up : number of neighbors is $vertex{$curvert}\n";
#	pick all the pairs containing picked vertex
	my %specedge=();
	foreach $m (keys %edge){
#uncomment for debug#		print "EDGES: $m  ----> @{$edge{$m}} \n";
		my @thisedge=@{$edge{$m}};
		$firstind=$[;
		$lastind=$#thisedge;
		$firstel=$thisedge[$firstind];
		$lastel=$thisedge[$lastind];
#uncomment for debug#		print "\$firstind=$firstind \$firstel=$firstel \$lastind=$lastind \$lastel=$lastel\n";
		# memorize the edge if it has $curvert at its ends
		if($firstel == $curvert && $lastel != $curvert){
			$specedge{$m}[0]='F';
			push @{$specedge{$m}},@thisedge;
		}
		if($lastel==$curvert && $firstel != $curvert){
			$specedge{$m}[0]='L';
			push @{$specedge{$m}},@thisedge;
		}

		if($lastel==$curvert && $firstel == $curvert){
			$specedge{$m}[0]='B';
			push @{$specedge{$m}},@thisedge;
		}
	}
#uncomment for debug#	foreach $m (keys %specedge){
#uncomment for debug#		print "SPECEDGES: $m  ----> @{$specedge{$m}} \n";
#uncomment for debug#	}
#	pick a number for a label; since we don't want it to accidentally repeat the existing labels, we choose
#       max of existing labels
	$maxlabel = -1;
	foreach $m (keys %edge){ $maxlabel=$m if $m > $maxlabel; }

#	join the two edges pairwise
	my @newedg=();
EDGE1:	foreach $m (keys %specedge){
		my @edge1=@{$specedge{$m}};
		next EDGE1 if $edge1[0] =~ 'B';
EDGE2:		foreach $n (keys %specedge){
		        my @edge2=@{$specedge{$n}};
			next EDGE2 if $edge2[0] =~ 'B';
			unless($m>=$n){
				@tmp1=@edge1;
				shift @tmp1;
				@tmp2=@edge2;
				shift @tmp2;
#uncomment for debug#				print "PAIRED SPECEDGES: $m  ----> @edge1 \n";
#uncomment for debug#				print "PAIRED SPECEDGES: $n  ----> @edge2 \n";
#uncomment for debug#				print "end of pair\n";
				# make the order of the first edge so that $curvert is last
				if($edge1[0] =~ 'F'){
					@tmp1= reverse @tmp1;
				}
				# make the order of the second edge so that $curvert is first
				if($edge2[0] =~ 'L'){
					@tmp2= reverse @tmp2;
				}
#uncomment for debug#				print "tmp1  ----> @tmp1 \n";
#uncomment for debug#				print "tmp2  ----> @tmp2 \n";
				shift @tmp2;
#uncomment for debug#				print "after tmp1  ----> @tmp1 \n";
#uncomment for debug#				print "after tmp2  ----> @tmp2 \n";
				# concatenate the two arrays together
				push(@tmp1,@tmp2);
#uncomment for debug#				print "joint tmp1  ----> @tmp1 \n";
				push @newedg, [ @tmp1 ];
			}
		}
	}
	for $m (0 .. $#newedg){
		$maxlabel+= 1;
        	for $j ( 0 .. $#{$newedg[$m]} ) {
            		#print "elementt $m $j is $newedg[$m][$j]\n";
			push @{$edge{$maxlabel}},$newedg[$m][$j];
        	}
	}
#uncomment for debug#	foreach $m (keys %edge){
#uncomment for debug#		print "NEW EDGES BEFORE REMOVING OLD ONES: $m  ----> @{$edge{$m}} \n";
#uncomment for debug#	}
CLEAN:	foreach $m (keys %edge){
		my @thisedge=@{$edge{$m}};
		$firstel=$thisedge[$[];
		$lastel=$thisedge[$#thisedge];
		if($firstel==$lastel){
			$ii=$#thisedge-1;
		}else{
			$ii=$#thisedge;
		}
		# remove edges which have repetitive elements along the path
		for $i (0 .. $ii){
			for $j (0 .. $i-1){
				if($thisedge[$j]==$thisedge[$i]){
					delete $edge{$m};
#uncomment for debug#					print "impossible ring path!: $thisedge[$i] is repetitive\n";
					next CLEAN;
				}
			}
		}
		if($firstel == $curvert || $lastel == $curvert){
			if($firstel == $lastel){
#uncomment for debug#				print "The ring is found!!! [ @thisedge ] \n";
				push @ring, [ @thisedge ];
			}
			delete $edge{$m};
			next CLEAN;
		}
	}
#uncomment for debug#	foreach $m (keys %edge){
#uncomment for debug#		print "NEW EDGES AFTER REMOVING THE OLD ONES: $m  ----> @{$edge{$m}} \n";
#uncomment for debug#	}
#	remove vertex
	delete $vertex{$curvert};
#       recalculate connectivities
#       zero out connectivities first
	foreach $m (keys %vertex){ 
		$vertex{$m}=0;
	}
	foreach $i (keys %vertex){ 
		foreach $m (keys %edge){
			my @thisedge=@{$edge{$m}};
			$firstel=$thisedge[$[];
			$lastel=$thisedge[$#thisedge];
			$vertex{$i}+=1 if ($firstel==$i || $lastel == $i);
		}
#uncomment for debug#			print "NEW VERTEX: $i $vertex{$i}\n";
	}

	} # end while

	foreach $m (0 .. $#ring){
#uncomment for debug#		print "Rings detected: @{$ring[$m]} \n";
		for $i ( 0 .. $#{$ring[$m]}-1 ){
			$tmp[0]=$ring[$m][$i];
			$tmp[1]=$ring[$m][$i+1];
			push @ringbonds, [ @ tmp ];
		}
	}
#uncomment for debug#	foreach $m (0 .. $#ringbonds){
#uncomment for debug#		print "RING BONDS: $m --> $ringbonds[$m][0] $ringbonds[$m][1]\n";
#uncomment for debug#	}

        my @keeplist=();
LABL:	for $k (0 .. $#m_bonds){ 
		$k1=$m_bonds[$k]->[1];
		$k2=$m_bonds[$k]->[2];
#uncomment for debug#				print "Analyze bond $k between $k1 and $k2\n";
		$t1=$types{$k1};
		$t2=$types{$k2};
#       check for bonds which are inside the rings 
		foreach $m (0 .. $#ringbonds){
			if(($k1==$ringbonds[$m][0] && $k2==$ringbonds[$m][1]) || ($k1==$ringbonds[$m][1] && $k2==$ringbonds[$m][0])){
#uncomment for debug#				print "Non-rotatable bond between $k1:$t1 and $k2:$t2 : in the ring\n";
				next LABL;
			}
		}
#       check for bonds which are connected on one end only
		$k1num= @{$neighb{$k1}};
		$k2num= @{$neighb{$k2}};
		if ($k1num==1||$k2num==1){
#uncomment for debug#			print "Non-rotatable bond between $k1:$t1 and $k2:$t2 : at the end\n";
			next LABL;
		}

#       check for imines 
		if($t1 =~ /N\.2/ && $t2 =~ /N\.2/ ){
#uncomment for debug#			print "Non-rotatable bond between $k1:$t1 and $k2:$t2 -- imine\n";
			next LABL;
		}
#       check for double bonds based on bond types
		for $m (0 .. $#m_bonds){ 
			$m1=$m_bonds[$m]->[1];
			$m2=$m_bonds[$m]->[2];
			$m3=$m_bonds[$m]->[3];
			if(($k1==$m1&&$k2==$m2)||($k1==$m2&&$k2==$m1)){
				if($m3 ne '1'){
#uncomment for debug#			                print "Non-rotatable bond between $k1:$t1 and $k2:$t2 -- non-single bond (bond type $m3)\n";
			                next LABL;
				}
			}
		}
#       check for double bonds,based on the atom types : for the sake of safety
		if($t1 =~ /C\.2/ && $t2 =~ /C\.2/ ){
#uncomment for debug#			print "Non-rotatable bond between $k1:$t1 and $k2:$t2 -- double bond\n";
			next LABL;
		}
#       check for triple bonds
		if($t1 =~ /C\.1/ && $t2 =~ /C\.1/ ){
#uncomment for debug#			print "Non-rotatable bond between $k1:$t1 and $k2:$t2 -- triple bond\n";
			next LABL;
		}
#       check for amides
		if( ($t1 =~ /C\.2/ && $t2 =~ /N\.2|N\.am|N\.2|N\.pl3/) || ($t2 =~ /C\.2/ && $t1 =~ /N\.2|N\.am|N\.2|N\.pl3/) ) {
#uncomment for debug#			print "Non-rotatable bond between $k1:$t1 and $k2:$t2 -- amide\n";
			next LABL;
		}
#       check for planar amines
		if( ($t1 =~ /C/ && $spectypes{$k2} =~ 'planar_nh2') ||  ($t2 =~ /C/ && $spectypes{$k1} =~ 'planar_nh2') ) {
#uncomment for debug#			print "Non-rotatable bond between $k1:$t1 and $k2:$t2 -- planar amine\n";
			next LABL;
		}
#       check for methyls
		if( ($t1 =~ /.*/ && $spectypes{$k2} =~ 'methyl') ||  ($t2 =~ /.*/ && $spectypes{$k1} =~ 'methyl') ) {
#uncomment for debug#			print "Non-rotatable bond between $k1:$t1 and $k2:$t2 -- methyl\n";
			next LABL;
		}
#       check for aromatic amidines
		if( ($t1 =~ /C\.ar/ && $spectypes{$k2} =~ 'amidine') ||  ($t2 =~ /C\.ar/ && $spectypes{$k1} =~ 'amidine') ) {
#uncomment for debug#			print "Non-rotatable bond between $k1:$t1 and $k2:$t2 -- aromatic amidine\n";
			next LABL;
		}
#       check for aromatic sulfonate
		if( ($t1 =~ /C\.ar/ && $spectypes{$k2} =~ 'sulfonate') ||  ($t2 =~ /C\.ar/ && $spectypes{$k1} =~ 'sulfonate') ) {
#uncomment for debug#			print "Non-rotatable bond between $k1:$t1 and $k2:$t2 -- aromatic amidine\n";
			next LABL;
		}
#       check for aromatic nitro/nitroso
		if( ($t1 =~ /C\.ar/ && $spectypes{$k2} =~ 'nitro') ||  ($t2 =~ /C\.ar/ && $spectypes{$k1} =~ 'nitro') ) {
#uncomment for debug#			print "Non-rotatable bond between $k1:$t1 and $k2:$t2 -- aromatic nitro/nitroso\n";
			next LABL;
		}
#       check for aromatic carboxylic acid. The big question remains if C.ar-COOH is rotatable - V.K.
		if( ($t1 =~ /C\.ar/ && $spectypes{$k2} =~ 'carboxy') ||  ($t2 =~ /C\.ar/ && $spectypes{$k1} =~ 'carboxy') ) {
		#uncomment for debug#	print "Non-rotatable bond between $k1:$t1 and $k2:$t2 -- aromatic carboxylic acid\n";
			next LABL;
		}
#	the remaining bonds should be considered rotatable 
		push @keeplist,$k;
	}
#uncomment for debug#	print "keeplist [@keeplist]\n";

	my @rotlist=();
	foreach $i (@keeplist) {
		$tmp[0]=$m_bonds[$i][1];
		$tmp[1]=$m_bonds[$i][2];
		push @rotlist, [@tmp] ;
	}    
	for $i (0 .. $#rotlist){
		print "Rotatable bond: $rotlist[$i][0] --- $rotlist[$i][1]\n";
	}
	$n_rbonds=@rotlist;
	print STDERR "number of rotatable bonds $n_rbonds\n";
	#my @metlist=();
	#foreach $i (keys %spectypes){
	#	if($spectypes{$i} =~ 'nos_methyl') {
	#		push(@metlist,$i);
	#		print "spectype for $i $spectypes{$i}\n";}
        #        }

        open(AUXFILE,">$file.aux") or die "Error while opening $file.aux $!\n";;
        print "Information about rotatable bonds will be output to $file.aux.\n";
	print AUXFILE "NROTBONDS  $n_rbonds\n";
	if($n_rbonds>0){
		for $i (0 .. $#rotlist){ print AUXFILE " $rotlist[$i][0] $rotlist[$i][1]\n"; }
	}
	close(AUXFILE);

	#$nummet = @metlist;
	#if($nummet > 0){
        #        print "Number of special methyls found: $nummet\n";
        #}

	close(MOL2FILE);
}
