// GMMDREC.cpp : Defines the entry point for the console application.
//
/********************************************************************/
/*  GMM-based recognizer                                            */
/*  Diagonal covariance matrices                                    */
/********************************************************************/
#include"stdafx.h"
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#define MAX_FEAT_LENGTH 40000				// Maximal Length of Feature Vector Sequence
#define MAX_PATTERN_NUM 500			        // Maximal Number of Classes to be identified
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
int MixtureNum,PatternNum;			
int FeatLen;						// Length of the Feature Vector Sequence (counted automatically )
int	ClipLen;
GMM *Model;        

int input(char *);
int pick_maximum(double *,int);
int cal_likelihood(int *,int,int);
double GaussianPDF(double *,double *,double *,double *,int);
void initialization(char *);
void Free(void);
void Exit(void);



void main(int argc,char *argv[])
{
	FILE	*fp1,*fp2;
	char	filename[200],FileAttribute[30];
	int		total,correct,ClipNum,i;
	int		result[500];

	if(argc<5){
		printf("Usage: GMMDREC TestCorpusFile GMMListFile PatternNum ResultFile ClipLen(#frames;'0':entire) ClipNum(1:only_one_clip)\n");
		exit(0);
	}

	PatternNum=atoi(argv[3]);
	initialization(argv[2]);
	if(argc>=6)	ClipLen=atoi(argv[5]);
	else	ClipLen=0;
	if(argc==7)	ClipNum=1;

	total=0;
	correct=0;
	fp2=fopen(argv[4],"w");
	if( (fp1=fopen(argv[1],"r"))!=NULL){
		while(fscanf(fp1,"%s %s\n",filename,FileAttribute)!=EOF){
			if((input(filename))!=-1){
				ClipNum=cal_likelihood(result,ClipLen,ClipNum);
				for(i=0;i<ClipNum;i++){
					total++;
					printf("%s %s  %s\n",filename,FileAttribute,ModelAttrib[result[i]]);
					fprintf(fp2,"%s %s  %s\n",filename,FileAttribute,ModelAttrib[result[i]]);
					if(strcmp(FileAttribute,ModelAttrib[result[i]])==0) correct++;
				}
			}
		}
	}
	printf("total=%d\n",total);
	printf("correct=%d\n",correct);
	printf("correct rate=%lf %%",(double)(correct)/(double)(total)*100.0);
	fprintf(fp2,"total=%d\n",total);
	fprintf(fp2,"correct=%d\n",correct);
	fprintf(fp2,"correct rate=%lf %%",(double)(correct)/(double)(total)*100.0);

	fclose(fp1);
	fclose(fp2);
	Free();
}



int input(char *filename)
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
	return 1;
}


int cal_likelihood(int *tend,int Len,int InNum)
{
	int		i,j,t,k,CNum;
	double	LogP[5000];	//max. mixture no. = 5000
	double	MaxLogP,temp;

	if(Len==0){
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

		tend[0]=pick_maximum(likelihood,PatternNum);
		return 1;
	}else{
		CNum=0;
		for(t=0;t<FeatLen;t +=Len){
			if((t+Len)>FeatLen) break;
			for(j=0;j<PatternNum;j++){
				likelihood[j]=0.0;
				for(k=0;k<Len;k++){
					MaxLogP=-100000000000000.0;
					for(i=0;i<MixtureNum;i++){
						LogP[i]=log(Model[j].GMMWeight[i])+GaussianPDF(Feat[t+k],Model[j].GMMMean[i],Model[j].GMMVar[i],Model[j].GMMVarconst[i],FeatDim);
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
//			printf("%d %.3f\n",j,likelihood[j]);
			tend[CNum]=pick_maximum(likelihood,PatternNum);

			CNum++;
			if(InNum==1)	break;
		}
		return CNum;
	}

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