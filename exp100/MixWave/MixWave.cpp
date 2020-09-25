// MixWave.cpp : Defines the entry point for the console application.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void main(short argc,char *argv[])
{
	double		energy1,energy2,EnergyRatio;
	short		data1,data2,y;
	FILE		*fp1,*fp2,*fp3;
	char		header1[44],header2[44];
	int			length1,length2,length,i;

	if(argc!=4){
		printf(" Usage : MixWave <Input WAV File 1> <Input WAV File 2> <Output WAV File>\n") ;
		exit(0) ;
	}

	if((fp1=fopen(argv[1],"rb"))==NULL){
		printf("%s error\n",argv[1]);
		exit(0);
	}
	if((fp2=fopen(argv[2],"rb"))==NULL){
		printf("%s error\n",argv[2]);
		exit(0);
	}
	if((fp3=fopen(argv[3],"wb"))==NULL){
		printf("disk space insufficient!\n");
		exit(0);
	}
	fread(header1,sizeof(char),44,fp1);
	fread(header2,sizeof(char),44,fp2);

/*
	fseek(fp1,0,SEEK_END);
	length1=(ftell(fp1)-44)/sizeof(short);
	fseek(fp1,44,0);
*/

	length1=0;
	energy1=0.0;
	while( fread(&data1,sizeof(short),1,fp1) == 1 ){  
		energy1 +=((double)data1)*((double)data1);
		length1++;
	}
	energy1 /=(double)length1;

	length2=0;
	energy2=0.0;
	while( fread(&data2,sizeof(short),1,fp2) == 1 ){    
		energy2 +=((double)data2)*((double)data2);
		length2++;
		if(length2==length1) break;
	}
	energy2 /=(double)length2;

	EnergyRatio=sqrt(energy1/energy2);

	if(length1<=length2){
		length=length1;
		fwrite(header1,sizeof(char),44,fp3);
	}else{
		length=length2;
		fwrite(header2,sizeof(char),44,fp3);
	}

	fseek(fp1,44,0);
	fseek(fp2,44,0);

	for(i=0;i<length;i++){
		fread(&data1,sizeof(short),1,fp1);    
		fread(&data2,sizeof(short),1,fp2);
		y=(short)(0.6/EnergyRatio*(double)(data1)+0.4*EnergyRatio*(double)(data2));
		fwrite(&y,sizeof(short),1,fp3);
	}

	fclose(fp1);
	fclose(fp2);
	fclose(fp3);
}

