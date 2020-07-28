#!/usr/bin/perl

#program to list the atoms that are dispalced from its position due to the 
#rotatable bonds and duplicates are removed.

if (@ARGV < 3) {
	die "Usage:$0 sdffile(s) auxfile(s) outputfile(s)file with .rot ext\n";
}

open (fsdf,$ARGV[0]) or die "sdf file not found";
@a_fsdf=<fsdf>;
#print @a_fsdf;

open (f,$ARGV[1]) or die "aux file not found";

open(fout,">".$ARGV[2]);

@a=<f>;
print (fout @a);
print (fout $a[0]);
foreach $line(@a)
{
	if ($line=~/^\s\d+\s(\d+)/) 
	{
		$rot=$1;
		@finalarray=sort(@finalarray);

		@finalarray=remove_duplicates(@finalarray);
		print (fout "Num_array_elment ".scalar(@finalarray)."\n");		

		print (fout "@finalarray\n");

		@finalarray=();
		foreach $sdf_line (@a_fsdf)
		{
		  #print $sdf_line;
		  if ($sdf_line=~/^\s+(\d+)\s+(\d+)\s+\d\s+\d$/)
		  {
			#print "$1,$2\n";
			$one=$1;$two=$2;	
			if ($1==$rot) 
			{
			  #print "$two-- \n"; 
			  push(@finalarray,$two);
			  Recursion1($two);
			}
		  }
		}
}
}
@finalarray=sort(@finalarray);
#print (fout "Num_array_elment ".scalar(@finalarray)."\n");		

		@finalarray=remove_duplicates(@finalarray);
		print (fout "Num_array_elment ".scalar(@finalarray)."\n");		

		print (fout "@finalarray\n");

sub Recursion1
{
my $var=@_[0]; 
foreach $sdf_line (@a_fsdf)
		{
		  #print $sdf_line;
		  if ($sdf_line=~/^\s+(\d+)\s+(\d+)\s+\d\s+\d$/)
		  {
			#print $1,$2;
			$one=$1;$two=$2;	
			if ($1==$var) { #print "$two "; 
 					push(@finalarray,$two);
					Recursion1($two);
				      }
		  }
		}
}


sub remove_duplicates
{
my (%seen, $finalarray);
($finalarray)=@_;
if(@$finalarray){return grep{$seen{$_} ++} @$finalarray;}
else{return grep{ ! $seen{$_} ++} @_;}
}



