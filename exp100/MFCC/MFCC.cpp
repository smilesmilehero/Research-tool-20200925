#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "wcomplex.h"

#define MAX_FILTERBANK_NO 22
#define	MAX_DIMENSION	100

short	FrameLength,FrameShift;
short	SamplingRate,Dimension,FilterBankNo;
float	f_begin[MAX_FILTERBANK_NO],f_end[MAX_FILTERBANK_NO],f_mid[MAX_FILTERBANK_NO];

void windowing(float *,float *,float *);
void Cal_frequency_tick(short,short,float *);
void Cal_triangular_bank_table(void);
void Cal_cep(float *,float *);
short filter_bank(float *,float *,float *,short);


void main(short argc,char *argv[])
{
	complex		*spec;
	float		y,tmp;
	float		*Data,*WindowedData,*hamming;
	float		*freq,*Magnitude;
	float		LogEnergy[MAX_FILTERBANK_NO],cep[MAX_DIMENSION];
	short		x,count;
	int			FeatLen=0;
	char		tempfile[200];
	FILE		*fp1,*fp2,*fp3;

	int			SampleStep=100000;				// dummy, for HTK format
  	short		SampleSize,SampleType=6;		// dummy, for HTK format

	if(argc!=5){
		printf(" Usage : MFCC <Input WAV File> <Output MFCC File> <SMagnitudeling Rate:8000/11025/16000/22050/etc.> <MFCC Dimension>\n") ;
		exit(0) ;
	}

	if((fp1=fopen(argv[1],"rb"))==NULL){
		printf("%s not found\n",argv[1]);
		exit(0);
	}
	fp2=fopen(argv[2],"wb") ;
	fseek(fp1,44,0);
	
	sprintf(tempfile,"%s.tmp",argv[1]);
	if((fp3=fopen(tempfile,"wb"))==NULL){
		printf("disk space insufficient!\n");
		exit(0);
	}

	SamplingRate=atoi(argv[3]);
	Dimension=atoi(argv[4]);
	SampleSize=Dimension*sizeof(float);

	if(SamplingRate<=8000)		FrameLength=256;
	if((SamplingRate>8000)&&(SamplingRate<=22050))	FrameLength=512;
	if(SamplingRate>22050)	FrameLength=1024;
	FrameShift=FrameLength/2;

	fwrite(&FeatLen,sizeof(int),1,fp2);
	fwrite(&SampleStep,sizeof(int),1,fp2);
	fwrite(&SampleSize,sizeof(short),1,fp2);
	fwrite(&SampleType,sizeof(short),1,fp2);

	Data=(float *)malloc(FrameLength*sizeof(float));
	WindowedData=(float *)malloc(FrameLength*sizeof(float));
	hamming=(float *)malloc(FrameLength*sizeof(float));
	freq=(float *)malloc(FrameLength*sizeof(float));
	Magnitude=(float *)malloc(FrameLength*sizeof(float));
	spec=(complex *)malloc(FrameLength*sizeof(complex));

	/******************* preemphasis ***************/
	/*   y[n]=x[n]-0.95*x[n-1]           ***********/
	/***********************************************/
	count=0;
	while( fread(&x,sizeof(short),1,fp1) == 1 ){ 
		if(count==0) y=(float)x;
		else y=(float)x-(float)0.97*tmp;
		tmp=(float)x;
		fwrite(&y,sizeof(float),1,fp3);
		count++;
	}
	fclose(fp3);

	/* hamming window table */
	for(count=0;count<FrameLength;count++)        
		hamming[count] = (float)(0.54 - 0.46 * cos( 2.0*3.141592*count / (FrameLength-1.0) ) );
	
	Cal_frequency_tick(SamplingRate,FrameLength,freq);
	
	Cal_triangular_bank_table();

	if(Dimension>FilterBankNo)	Dimension=FilterBankNo;

//	printf("FilterBankNo=%d, Dimension=%d\n",FilterBankNo,Dimension);

	fp3=fopen(tempfile,"rb");
	FeatLen=0;
	while( (short)fread(Data,sizeof(float),FrameLength,fp3) == FrameLength ){

		windowing(Data,WindowedData,hamming) ;

		for(count=0;count<FrameLength;count++){
			spec[count].re=WindowedData[count];
			spec[count].im=0.0;
		}
      
		fft(spec,FrameLength);

		for(count=0;count<FrameLength;count++)
			Magnitude[count]=(float)(sqrt(spec[count].re*spec[count].re+spec[count].im*spec[count].im));

		if( (filter_bank(Magnitude,LogEnergy,freq,FrameLength))==1 ){
			Cal_cep(LogEnergy,cep);
			fwrite(cep,sizeof(float),Dimension,fp2);
		}
		FeatLen++;

		fseek(fp3,FrameShift*FeatLen*sizeof(float),SEEK_SET);
	}

	fseek(fp2,0,SEEK_SET);
	fwrite(&FeatLen,sizeof(int),1,fp2);
	fwrite(&SampleStep,sizeof(int),1,fp2);
	SampleSize=Dimension*sizeof(float);
	fwrite(&SampleSize,sizeof(short),1,fp2);


	free(Data);
	free(WindowedData);
	free(hamming);
	free(freq);
	free(Magnitude);
	free(spec);

	fclose(fp1);
	fclose(fp2);
	fclose(fp3);
	remove(tempfile);
}

void windowing(float *d,float *dw,float *ham)
{
   short count ;

   for (count=0;count<FrameLength;count++)
       *(dw+count) = *(d+count) * *(ham+count) ;
}


void Cal_frequency_tick(short samplingRate,short size,float *freq)
{
  short k;

  for(k=0;k<size/2;k++)
    freq[k]=(float)k/(float)size*(float)samplingRate;
  for(k=size/2;k<size;k++)
    freq[k]=(float)(size-k)/(float)size*(float)samplingRate;

}


float mel_to_freq(float mel)
{
  return (float)(700.0*(pow(10.0,mel/2595.0)-1.0));
}


void Cal_triangular_bank_table(void)
{
  short i;            
  float	MelInterval=2590.0*log10(1+11025.0/700.0)/(MAX_FILTERBANK_NO);  // ¥Ñ FILTERBANK_NO ¨M©w MelInterval
  float m_begin[MAX_FILTERBANK_NO],m_end[MAX_FILTERBANK_NO],m_mid[MAX_FILTERBANK_NO];

  for(i=0;i<MAX_FILTERBANK_NO;i++){
    m_begin[i]=(float)i*MelInterval;
    m_end[i]=m_begin[i]+MelInterval*(float)2.0;
    m_mid[i]=m_begin[i]+MelInterval;
  }
  for(i=0;i<MAX_FILTERBANK_NO;i++){
    f_begin[i]=mel_to_freq(m_begin[i]);
    f_end[i]=mel_to_freq(m_end[i]);
    f_mid[i]=mel_to_freq(m_mid[i]);

	if(f_end[i]>=(float)SamplingRate/2)	break;
  }
  FilterBankNo=i;
 // printf("%hd\n",i);
}

float triangular(float begin,float end,float mid,float x)
{
	if( (x<mid)&&(x>begin) )    return (x-begin)/(mid-begin);
	if( (x>mid)&&(x<end) )    return (x-end)/(mid-end);
	if(x==mid)    return 1.0;
    
	return 0.0;
}


short filter_bank(float *magnitude,float *log_energy,float *freq,short size)
{
	short	i,j;
	float	response;
	float	df;

	df=freq[1]-freq[0];

	for(i=0;i<FilterBankNo;i++){
		log_energy[i]=0.0;
		for(j=0;j<size;j++){
			response=magnitude[j]*triangular(f_begin[i],f_end[i],f_mid[i],freq[j]);
			log_energy[i] +=response*response*df;
		}

		if(log_energy[i]==0.0) return 0;
		else log_energy[i]=(float)log(log_energy[i]);
	}
	
	return 1;
}

void Cal_cep(float *E,float *C)
{
  short k,m;

  for(m=0;m<Dimension;m++){
    C[m]=0.0;
    for(k=0;k<FilterBankNo;k++)
      C[m] +=E[k]*(float)(cos( (float)(m+1.0)*(float)(k+0.5)*3.14159/(float)FilterBankNo ));

//	printf("%d %f\n",m,C[m]);
  }
 // for(m=0;m<Dimension;m++){
//	  if(m>8) {
//	  C[m]=0.0;
//	  }
 // else
//	  continue;
 // }
}

