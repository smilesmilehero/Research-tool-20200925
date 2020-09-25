/********************************************************************/
/*  GMM-based recognizer                                            */
/*  Diagonal covariance matrices                                    */
/*  Training + testing                                              */
/********************************************************************/
#include"stdafx.h"
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include <process.h>
#define MAX_FEAT_LENGTH 40000				// Maximal Length of Feature Vector Sequence
#define MAX_PATTERN_NUM 2000			        // Maximal Number of Classes to be identified
typedef struct{
	double	**GMMMean;
	double	**GMMVar;
	double	**GMMVarconst;
	double	*GMMWeight;
}GMM;

double **Feat;									// Feature Vector Feat[Frame_no][dimension]
double likelihood[MAX_PATTERN_NUM];
char ModelAttrib[MAX_PATTERN_NUM][30];			
int FeatDim;									// Dimensionality of Feature Vector
int MixtureNum=16,PatternNum,SpeakerNum;			
int FeatLen;						// Length of the Feature Vector Sequence (counted automatically )
GMM *Model;        

void input(char *);
int pick_maximum(double *,int);
int cal_likelihood(void);
double GaussianPDF(double *,double *,double *,double *,int);
void initialization(char *);
void Free(void);
void Exit(void);



void main(int argc,char *argv[])
{
	FILE	*fp;
	int		SimultaneousSpeakerNo=1,IfTrain=1,result,i,j,k,m;
	char	string1[30],string2[30],string3[30],string4[30],string5[30];

	if(argc!=4){
		printf("Usage: GMMSPID Simultaneous-Speaker-Number(default=1) IfTrain(default=1) MixtureNum(default=16)\n");
		exit(0);
	}

	SimultaneousSpeakerNo=atoi(argv[1]);
	IfTrain=atoi(argv[2]);
	MixtureNum=atoi(argv[3]);

	printf("Enter the number of speakers (from 2 to 10):  ");
	scanf("%d",&SpeakerNum);
	printf("%d Speakers\n",SpeakerNum);

	if(IfTrain==1){
		printf("Please wait for training...\n");

		/******* single-speaker GMM training *******/
		fp=fopen("train-1speaker.crp","w");
		for(i=0;i<SpeakerNum;i++){
			sprintf(string1,"wav\\%02d-01.wav",i);
			sprintf(string2,"mfc\\%02d-01.mfc",i);
			_spawnl(_P_WAIT,"mfcc","mfcc",string1,string2,"8000","15",NULL);
			fprintf(fp,"%s %02d\n",string2,i);
		}
		fclose(fp);
	
		fp=fopen("GMM-1Speaker.lst","w");
		for(i=0;i<SpeakerNum;i++){
			sprintf(string1,"ref\\GMM-%02d",i);
			sprintf(string2,"%02d",i);
			_spawnl(_P_WAIT,"gmmd","gmmd","train-1speaker.crp",string1,argv[3],string2,"20",NULL);
			fprintf(fp,"%s %s\n",string1,string2);
		}
		fclose(fp);
	
	
		if(SimultaneousSpeakerNo>=2){	/******* two-speaker GMM training ******/

			for(i=0;i<SpeakerNum;i++){			// mixing waves 
				for(j=0;j<SpeakerNum;j++){
					if(i!=j){
						sprintf(string1,"wav\\%02d-01.wav",i);
						sprintf(string2,"wav\\%02d-01.wav",j);
						sprintf(string3,"wav\\%02d-%02d-01.wav",i,j);
						_spawnl(_P_WAIT,"MixWave","MixWave",string1,string2,string3,NULL);
					}
				}
			}

			fp=fopen("train-2speakers.crp","w");
			for(i=0;i<SpeakerNum;i++){
				for(j=0;j<SpeakerNum;j++){
					if(i!=j){
						sprintf(string1,"wav\\%02d-%02d-01.wav",i,j);
						sprintf(string2,"mfc\\%02d-%02d-01.mfc",i,j);
						_spawnl(_P_WAIT,"mfcc","mfcc",string1,string2,"8000","15",NULL);
						fprintf(fp,"%s %02d-%02d\n",string2,i,j);
					}
				}
			}
			fclose(fp);

			fp=fopen("GMM-2Speakers.lst","w");
			for(i=0;i<SpeakerNum;i++){
				for(j=0;j<SpeakerNum;j++){
					if(i!=j){
						sprintf(string1,"ref\\GMM-%02d-%02d",i,j);
						sprintf(string2,"%02d-%02d",i,j);
						_spawnl(_P_WAIT,"gmmd","gmmd","train-2speakers.crp",string1,argv[3],string2,"20",NULL);
					}else{
						sprintf(string1,"ref\\GMM-%02d",i);
						sprintf(string2,"%02d",i);
					}
					fprintf(fp,"%s %s\n",string1,string2);
				}
			}
			fclose(fp);
		}


		if(SimultaneousSpeakerNo>=3){	/******* three-speaker GMM training ******/

			for(i=0;i<SpeakerNum;i++){			// mixing waves 
				for(j=0;j<SpeakerNum;j++){
					for(k=0;k<SpeakerNum;k++){
						if((i!=j)&&(k!=i)&&(k!=j)){
							sprintf(string3,"wav\\%02d-%02d-01.wav",i,j);
							sprintf(string4,"wav\\%02d-01.wav",k);
							sprintf(string5,"wav\\%02d-%02d-%02d-01.wav",i,j,k);
							_spawnl(_P_WAIT,"MixWave","MixWave",string3,string4,string5,NULL);
						}
					}
				}
			}

			fp=fopen("train-3speakers.crp","w");
			for(i=0;i<SpeakerNum;i++){
				for(j=0;j<SpeakerNum;j++){
					for(k=0;k<SpeakerNum;k++){
						if((i!=j)&&(k!=i)&&(k!=j)){
							sprintf(string1,"wav\\%02d-%02d-%02d-01.wav",i,j,k);
							sprintf(string2,"mfc\\%02d-%02d-%02d-01.mfc",i,j,k);
							_spawnl(_P_WAIT,"mfcc","mfcc",string1,string2,"8000","15",NULL);
							fprintf(fp,"%s %02d-%02d-%02d\n",string2,i,j,k);
						}
					}
				}
			}
			fclose(fp);

			fp=fopen("GMM-3Speakers.lst","w");
			for(i=0;i<SpeakerNum;i++){
				for(j=0;j<SpeakerNum;j++){
					for(k=0;k<SpeakerNum;k++){
						if((i!=j)&&(k!=i)&&(k!=j)){
							sprintf(string1,"ref\\GMM-%02d-%02d-%02d",i,j,k);
							sprintf(string2,"%02d-%02d-%02d",i,j,k);
							_spawnl(_P_WAIT,"gmmd","gmmd","train-3speakers.crp",string1,argv[3],string2,"20",NULL);
						}else if((i!=j)&&(k==j)){
							sprintf(string1,"ref\\GMM-%02d-%02d",i,j);
							sprintf(string2,"%02d-%02d",i,j);
						}else if((i==j)&&(k!=j)){
							sprintf(string1,"ref\\GMM-%02d-%02d",j,k);
							sprintf(string2,"%02d-%02d",j,k);
						}else if((i==j)&&(i==k)){
							sprintf(string1,"ref\\GMM-%02d",i);
							sprintf(string2,"%02d",i);
						}
						fprintf(fp,"%s %s\n",string1,string2);
					}
				}
			}	
			fclose(fp);
		}

		if(SimultaneousSpeakerNo>=4){	/******** four-speaker GMM training *******/

			for(i=0;i<SpeakerNum;i++){			// mixing waves 
				for(j=0;j<SpeakerNum;j++){
					for(k=0;k<SpeakerNum;k++){
						for(m=0;m<SpeakerNum;m++){
							if((i!=j)&&(k!=i)&&(k!=j)&&(m!=i)&&(m!=j)&&(m!=k)){
								sprintf(string3,"wav\\%02d-%02d-%02d-01.wav",i,j,k);
								sprintf(string4,"wav\\%02d-01.wav",m);
								sprintf(string5,"wav\\%02d-%02d-%02d-%02d-01.wav",i,j,k,m);
								_spawnl(_P_WAIT,"MixWave","MixWave",string3,string4,string5,NULL);
							}
						}
					}
				}
			}

			fp=fopen("train-4speakers.crp","w");
			for(i=0;i<SpeakerNum;i++){
				for(j=0;j<SpeakerNum;j++){
					for(k=0;k<SpeakerNum;k++){
						for(m=0;m<SpeakerNum;m++){
							if((i!=j)&&(k!=i)&&(k!=j)&&(m!=i)&&(m!=j)&&(m!=k)){
								sprintf(string1,"wav\\%02d-%02d-%02d-%02d-01.wav",i,j,k,m);
								sprintf(string2,"mfc\\%02d-%02d-%02d-%02d-01.mfc",i,j,k,m);
								_spawnl(_P_WAIT,"mfcc","mfcc",string1,string2,"8000","15",NULL);
								fprintf(fp,"%s %02d-%02d-%02d-%02d\n",string2,i,j,k,m);
							}
						}
					}
				}
			}
			fclose(fp);

			fp=fopen("GMM-4Speakers.lst","w");
			for(i=0;i<SpeakerNum;i++){
				for(j=0;j<SpeakerNum;j++){
					for(k=0;k<SpeakerNum;k++){
						for(m=0;m<SpeakerNum;m++){
							if((i!=j)&&(k!=i)&&(k!=j)&&(m!=i)&&(m!=j)&&(m!=k)){
								sprintf(string1,"ref\\GMM-%02d-%02d-%02d-%02d",i,j,k,m);
								sprintf(string2,"%02d-%02d-%02d-%02d",i,j,k,m);
								_spawnl(_P_WAIT,"gmmd","gmmd","train-4speakers.crp",string1,argv[3],string2,"20",NULL);

							}else if((i==j)&&(k!=i)&&(m!=i)&&(m!=k)){
								sprintf(string1,"ref\\GMM-%02d-%02d-%02d",i,k,m);
								sprintf(string2,"%02d-%02d-%02d",i,k,m);
							}else if((i==k)&&(j!=i)&&(m!=i)&&(m!=j)){
								sprintf(string1,"ref\\GMM-%02d-%02d-%02d",i,j,m);
								sprintf(string2,"%02d-%02d-%02d",i,j,m);
							}else if((i==m)&&(j!=i)&&(k!=i)&&(j!=k)){
								sprintf(string1,"ref\\GMM-%02d-%02d-%02d",i,j,k);
								sprintf(string2,"%02d-%02d-%02d",i,j,k);
							}else if((j==k)&&(i!=j)&&(m!=j)&&(m!=i)){
								sprintf(string1,"ref\\GMM-%02d-%02d-%02d",i,j,m);
								sprintf(string2,"%02d-%02d-%02d",i,j,m);
							}else if((j==m)&&(i!=m)&&(k!=m)&&(i!=k)){
								sprintf(string1,"ref\\GMM-%02d-%02d-%02d",i,j,k);
								sprintf(string2,"%02d-%02d-%02d",i,j,k);
							}else if((k==m)&&(i!=k)&&(j!=k)&&(i!=j)){
								sprintf(string1,"ref\\GMM-%02d-%02d-%02d",i,j,k);
								sprintf(string2,"%02d-%02d-%02d",i,j,k);

							}else if((i==j)&&(i==k)&&(m!=i)){
								sprintf(string1,"ref\\GMM-%02d-%02d",i,m);
								sprintf(string2,"%02d-%02d",i,m);
							}else if((i==j)&&(i==m)&&(k!=i)){
								sprintf(string1,"ref\\GMM-%02d-%02d",i,k);
								sprintf(string2,"%02d-%02d",i,k);
							}else if((i==k)&&(i==m)&&(j!=i)){
								sprintf(string1,"ref\\GMM-%02d-%02d",i,j);
								sprintf(string2,"%02d-%02d",i,j);

							}else if((i==j)&&(i==k)&&(i==m)){
								sprintf(string1,"ref\\GMM-%02d",i);
								sprintf(string2,"%02d",i);
							}
							fprintf(fp,"%s %s\n",string1,string2);
						}
					}
				}
			}	
			fclose(fp);
		}
	}

	switch(SimultaneousSpeakerNo){
		case 1:
			PatternNum=SpeakerNum;
			initialization("GMM-1Speaker.lst");
			break;
		case 2:
			PatternNum=SpeakerNum*SpeakerNum;
			initialization("GMM-2Speakers.lst");
			break;
		case 3:
			PatternNum=SpeakerNum*SpeakerNum*SpeakerNum;
			initialization("GMM-3Speakers.lst");
			break;
		case 4:
			PatternNum=SpeakerNum*SpeakerNum*SpeakerNum*SpeakerNum;
			initialization("GMM-4Speakers.lst");
			break;

	}
	

	/***** Identification ***************/
	printf("Press <Enter> to perform identification, or press <q> to escape.\n");
	getc(stdin);
	while((getc(stdin))!='q'){
		_spawnl(_P_WAIT,"mfcc","mfcc","test.wav","test.mfc","8000","15",NULL);
	
		input("test.mfc");		
	
		result=cal_likelihood();

		printf("Identified Speaker: %s\n",ModelAttrib[result]);
	
		printf("Press <Enter> to perform identification, or press <q> to escape.\n");

	}


	Free();
}



void input(char *filename)
{
	int		i,t,SampleStep;
	FILE	*fpFeat;
  	short	SampleSize,SampleType;
	float	temp[100];

	if((fpFeat=fopen(filename,"rb"))==NULL){
		printf("%s opened error\n",filename);
		exit(0);
	}

	fread(&FeatLen,sizeof(int),1,fpFeat);
	fread(&SampleStep,sizeof(int),1,fpFeat);
	fread(&(SampleSize),sizeof(short),1,fpFeat);
	if(int(SampleSize/sizeof(float))!=FeatDim){
		printf("dimension mismatch between features and models\n");
		exit(0);
	}
	fread(&SampleType,sizeof(short),1,fpFeat);
	FeatDim=SampleSize/sizeof(float);
	for(t=0;t<FeatLen;t++){
		if(int(fread(temp,sizeof(float),FeatDim,fpFeat))==FeatDim){
			for(i=0;i<FeatDim;i++){
				Feat[t][i]=(double)temp[i];
			}
		}
	}
	fclose(fpFeat);

}


int cal_likelihood(void)
{
	int		i,j,t,tend;
	double	LogP[5000];	//max. mixture no. = 5000
	double	MaxLogP,temp;

	for(j=0;j<PatternNum;j++){	
			likelihood[j]=0.0;
			for(t=0;t<FeatLen;t++){

				MaxLogP=-100000000000000.0;
				for(i=0;i<MixtureNum;i++){
					LogP[i]=log(Model[j].GMMWeight[i])+GaussianPDF(Feat[t],Model[j].GMMMean[i],Model[j].GMMVar[i],Model[j].GMMVarconst[i],FeatDim);
					if(LogP[i]>MaxLogP)	MaxLogP=LogP[i];
				}
				temp=0.0;
				for(i=0;i<MixtureNum;i++)	
					if((LogP[i]-MaxLogP)>-744.0)	
						temp +=exp(LogP[i]-MaxLogP);

				if(temp!=0.0)  likelihood[j] +=log(temp);
				likelihood[j] +=MaxLogP;

			}
	}
	
	tend=pick_maximum(likelihood,PatternNum);

	return tend;

}


int pick_maximum(double *q,int numb)
{
  int i,j;
  double temp;

  temp=q[0];
  j=0;
  for(i=1;i<numb;i++){
	if(q[i]>=temp){
      j=i;
	  temp=q[i];
	}
  }

  return j;
}

double GaussianPDF(double *GX,double *MeAn,double *VaR,double *VaRcOnSt,int Dim)
{
  double p=0.0;
  int i;

  for(i=0;i<Dim;i++)
    p += VaRcOnSt[i]+((GX[i]-MeAn[i])*(GX[i]-MeAn[i])/(-2.0*VaR[i])) ;

  return p;
}


void initialization(char *filename)
{
	FILE *fp,*fp2;
	char GMMFile[200];
	int i,j,k,index;

	if((fp=fopen(filename,"r"))==NULL){
		printf("%s not found\n",filename);
		exit(0);
	}

	Model=(GMM *)malloc(PatternNum*sizeof(GMM));

	index=0;
	while((fscanf(fp,"%s %s\n",GMMFile,ModelAttrib[index]))!=EOF){
		if((fp2=fopen(GMMFile,"rb"))!=NULL){
			fread(&FeatDim,sizeof(int),1,fp2);
			fread(&MixtureNum,sizeof(int),1,fp2);
			Model[index].GMMMean=(double **)malloc(MixtureNum*sizeof(double *));
			Model[index].GMMVar=(double **)malloc(MixtureNum*sizeof(double *));
			Model[index].GMMVarconst=(double **)malloc(MixtureNum*sizeof(double *));
			Model[index].GMMWeight=(double *)malloc(MixtureNum*sizeof(double));
			for(i=0;i<MixtureNum;i++){
				Model[index].GMMMean[i]=(double *)malloc(FeatDim*sizeof(double));
				Model[index].GMMVar[i]=(double *)malloc(FeatDim*sizeof(double));
				Model[index].GMMVarconst[i]=(double *)malloc(FeatDim*sizeof(double));
			}
			for(i=0;i<MixtureNum;i++)
				fread(Model[index].GMMMean[i],sizeof(double),FeatDim,fp2);
			for(i=0;i<MixtureNum;i++)
				fread(Model[index].GMMVar[i],sizeof(double),FeatDim,fp2);
			fread(Model[index].GMMWeight,sizeof(double),MixtureNum,fp2);
			fclose(fp2);
			index++;
		}else{
			printf("%s not found\n",GMMFile);
			exit(0);
		}
	}
	fclose(fp);

	for(k=0;k<PatternNum;k++){
		for(i=0;i<MixtureNum;i++){
			for(j=0;j<FeatDim;j++){
				if(Model[k].GMMVar[i][j]==0.0){
					printf("covariance=0\n");
					exit(0);
				}
				Model[k].GMMVarconst[i][j]=log(1.0/sqrt(Model[k].GMMVar[i][j]*6.28318));
			}
		}
	}  

	if((Feat=(double **)malloc(MAX_FEAT_LENGTH*sizeof(double *)))==NULL) Exit();
	for(i=0;i<MAX_FEAT_LENGTH;i++)
		if((Feat[i]=(double *)malloc(FeatDim*sizeof(double)))==NULL) Exit();
}

void Free(void)
{
	int i,j;

	for(i=0;i<MAX_FEAT_LENGTH;i++)	free(Feat[i]);
	free(Feat);
	for(j=0;j<PatternNum;j++){
		for(i=0;i<MixtureNum;i++){
			free(Model[j].GMMMean[i]);
			free(Model[j].GMMVar[i]);
			free(Model[j].GMMVarconst[i]);
		}
		free(Model[j].GMMMean);
		free(Model[j].GMMVar);
		free(Model[j].GMMVarconst);
		free(Model[j].GMMWeight);
	}
	free(Model);
}

void Exit(void)
{
  printf("memory insufficient\n");
  exit(0);
}