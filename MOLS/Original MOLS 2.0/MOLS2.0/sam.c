#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

void sort( int [], int );
void sort( int a[], int elements )
{
   int i, j, temp;
   i = 0;
   while( i < (elements - 1) )
  {
   j = i + 1;
   while( j < elements )
  {
   if( a[i] > a[j] )
    {
    temp = a[i];
    a[i] = a[j];
    a[j] = temp;
    }
  j++;
  }
 i++;
  }
}

main(int argc, char *argv[])
{
char str[20];
int count=1,totno,nrot,i,j,k,tmp,loop;
int at1[100],at2[100];
int numb[1000];

FILE *fin,*fout;

	if( argc != 3 )
         {
         printf("Incorrect, Useage: <./a.out> <rot file> <out file> \n");
         exit(2);
         }

fin= fopen(argv[1],"r");
	if (fin==NULL){
	printf("File does not exists\n");
	exit(0);
	}

fout=fopen(argv[2],"w");

	fscanf(fin,"%s%d",str,&nrot);
	fprintf(fout,"%9s %3d\n",str,nrot);
	for (i=1;i<=nrot;i++)
	{
	fscanf(fin,"%d%d",&at1[i],&at2[i]);
//	fprintf(fout,"%d\t%d\n",at1,at2);
	}
	fscanf(fin,"%s%d",str,&nrot);
	fscanf(fin,"%s%d",str,&totno);
	fscanf(fin,"\n");
	for (i=1;i<=nrot;i++)
	{
	fprintf(fout,"%3d     %3d\n",at1[i],at2[i]);
	fscanf(fin,"%s%d",str,&totno);
	fprintf(fout,"%16s %3d\n",str,totno);

		for (j=0;j<totno;j++)
		{
		fscanf(fin,"%d",&numb[j]);
		}
	sort( numb,totno);
       for( loop = 0; loop<totno; loop++ )
       fprintf(fout,"%3d\n",numb[loop]);
	}


fclose(fin);
fclose(fout);


}




