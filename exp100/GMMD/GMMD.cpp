// GMMD.cpp : Defines the entry point for the console application.
//
/********************************************************************/
/* Generating a Gaussian Mixture Model (GMM)	         			*/
/* Diagonal-covariance matrices.									*/
/* GMM parameters are initialized via K-means Algorithm             */
/* and are reestimated via EM Algorithm                             */
/********************************************************************/
#include"stdafx.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#define MAX_TRA_FILE_NUM 500      // Maximal Number of Training Files
#define MAX_ITER_NUM 50           // Maximal Iteration No. in EM

float ***Feat;                    // Feature Vector Feat[File_no][Frame_no][Element_no]
double **GMMMean;                 // Mean[Mixture_no][Element_no];
double **GMMVar;                  // Variance[Mixture_no][Element_no];
double *GMMWeight;                // Weight[Mixture_no];
double **GMMVarconst;			  // denominator of Gaussian pdf
int MixtureNum;                   // given from argument
int FileNum;					  // automatically count
int FeatDim;				      // Dimensionality of Feature Vector
int IterationNum;				  // given from argument
int FeatLen[MAX_TRA_FILE_NUM],TotalLen;  // FeatLen[i] is the length in ith file, and TotalLen is the total length of all training files.

double GaussianPDF(float *,double *,double *,double *,int);
double distance(float *,double *,int);
void get_training_files(char *,char *);
void initialization(void);
void K_means(void);
void reestimation(char *);
void save_GMM(char *);
void Free(void);
void Exit(void);
int pick_minimum(double *,int);


void main(int argc,char *argv[])
{

  if(argc!=6){
	printf("Usage: GMMD CorpusFile SaveGMMFile MixtureNum Identity IterationNo\n");
	exit(0);
  }

  IterationNum=atoi(argv[5]);
  get_training_files(argv[1],argv[4]);
  MixtureNum=atoi(argv[3]);
  initialization();
  K_means();
  reestimation(argv[2]);
  Free();

}

void get_training_files(char *datfile,char *target)
{
	int		f,t,SampleStep;
	FILE	*fpDat,*fpFeat;
	char	filename[200],Identity[30];
  	short	SampleSize,SampleType;

	if((Feat=(float ***)malloc(MAX_TRA_FILE_NUM*sizeof(float **)))==NULL) Exit();

	f=0;
	TotalLen=0;
	if( (fpDat=fopen(datfile,"r"))!=NULL){
		while(fscanf(fpDat,"%s %s\n",filename,Identity)!=EOF){
			if(strcmp(Identity,target)==0){
				if((fpFeat=fopen(filename,"rb"))==NULL){
					printf("%s opened error\n",filename);
					exit(0);
				}
				fread(&(FeatLen[f]),sizeof(int),1,fpFeat);
				fread(&SampleStep,sizeof(int),1,fpFeat);
				fread(&(SampleSize),sizeof(short),1,fpFeat);
				fread(&SampleType,sizeof(short),1,fpFeat);
				FeatDim=SampleSize/sizeof(float);
//				printf("%s  %d\n",filename,FeatLen[f]);
				if((Feat[f]=(float **)malloc(FeatLen[f]*sizeof(float *)))==NULL)  Exit();
				for(t=0;t<FeatLen[f];t++){
					if((Feat[f][t]=(float *)malloc(FeatDim*sizeof(float)))==NULL)  Exit();
					fread(Feat[f][t],sizeof(float),FeatDim,fpFeat);
				}
				TotalLen=TotalLen+FeatLen[f];
				f++;
				fclose(fpFeat);
			}
		}
	}else{
		printf("%s was not found!",datfile);
		exit(0);
	}
	FileNum=f;
//	printf("number of training files=%d\n",FileNum);
	if(FileNum==0) exit(0);
//	printf("Speaker %s is trained OK\n",target);

	fclose(fpDat);
}




void initialization(void)
{
  int i;

  if((GMMMean=(double **)malloc(MixtureNum*sizeof(double *)))==NULL) Exit();
  if((GMMVar=(double **)malloc(MixtureNum*sizeof(double *)))==NULL) Exit();
  if((GMMWeight=(double *)malloc(MixtureNum*sizeof(double)))==NULL) Exit();
  if((GMMVarconst=(double **)malloc(MixtureNum*sizeof(double *)))==NULL) Exit();
  for(i=0;i<MixtureNum;i++){
    if((GMMMean[i]=(double *)malloc(FeatDim*sizeof(double)))==NULL)      Exit();
    if((GMMVar[i]=(double *)malloc(FeatDim*sizeof(double)))==NULL)       Exit();
    if((GMMVarconst[i]=(double *)malloc(FeatDim*sizeof(double)))==NULL)       Exit();
  }

}


void reestimation(char *savefile)
{
	int		i,j,k,t,f;
	double	***likeli,**likeli_sum;
	double	LogP[5000],MaxLogP;			// max. Mix. no. = 5000
	double	OverallScore[MAX_ITER_NUM],temp1[100],temp3[100],temp2,temp4;

	for(i=0;i<MixtureNum;i++)  GMMWeight[i]=1.0/(double)MixtureNum;
	if((likeli=(double ***)malloc(FileNum*sizeof(double **)))==NULL) Exit();
	for(f=0;f<FileNum;f++){
		if((likeli[f]=(double **)malloc(FeatLen[f]*sizeof(double *)))==NULL)  Exit();
		for(t=0;t<FeatLen[f];t++)
		if((likeli[f][t]=(double *)malloc(MixtureNum*sizeof(double)))==NULL)  Exit();
	}
	if((likeli_sum=(double **)malloc(FileNum*sizeof(double *)))==NULL) Exit();
	for(f=0;f<FileNum;f++)	if((likeli_sum[f]=(double *)malloc(FeatLen[f]*sizeof(double)))==NULL)  Exit();

	// maximum OverallScore search
	for(k=0;k<IterationNum;k++){

		OverallScore[k]=0.0;
		for(f=0;f<FileNum;f++){
			for(t=0;t<FeatLen[f];t++){
				
				MaxLogP=-10000000000000.0;
				for(j=0;j<MixtureNum;j++){
					LogP[j]=log(GMMWeight[j])+GaussianPDF(Feat[f][t],GMMMean[j],GMMVar[j],GMMVarconst[j],FeatDim);
					if(LogP[j]>MaxLogP)		MaxLogP=LogP[j];
				}

				likeli_sum[f][t]=0.0;
				for(j=0;j<MixtureNum;j++){
					if((LogP[j]-MaxLogP)>-744.0)	likeli[f][t][j]=exp(LogP[j]-MaxLogP);
					else	likeli[f][t][j]=0.0;
					likeli_sum[f][t] +=likeli[f][t][j];
				}
				if(likeli_sum[f][t]!=0.0)	OverallScore[k] +=log(likeli_sum[f][t]);
				OverallScore[k] +=MaxLogP;
			}
		}
//		printf("%d: %lf ",k,OverallScore[k]);

		for(i=0;i<MixtureNum;i++){
			for(j=0;j<FeatDim;j++)	temp1[j]=temp3[j]=0.0;
			temp2=0.0;
			for(f=0;f<FileNum;f++){
				for(t=0;t<FeatLen[f];t++){
					if(likeli_sum[f][t]!=0.0)	temp4=likeli[f][t][i]/likeli_sum[f][t];
					else	temp4=0.0;
					temp2 +=temp4;
					for(j=0;j<FeatDim;j++){
						temp1[j] +=temp4*(double)(Feat[f][t][j]);
						temp3[j] +=temp4*((double)(Feat[f][t][j])-GMMMean[i][j])*((double)(Feat[f][t][j])-GMMMean[i][j]);
					}
				}
			}
			for(j=0;j<FeatDim;j++){
				GMMMean[i][j]=temp1[j]/temp2;
				GMMVar[i][j]=temp3[j]/temp2;
	  			if(GMMVar[i][j]<0.00000001)		GMMVar[i][j]=0.00000001;
				GMMVarconst[i][j]=log(1.0/sqrt(GMMVar[i][j]*6.28318));
			}
			GMMWeight[i]=temp2/(double)(TotalLen);
		}

		if( (OverallScore[k]>=OverallScore[k-1])&&(OverallScore[k]-OverallScore[k-1]<0.01)) break;
		if( (k>10)&&(OverallScore[k]<=OverallScore[k-1]) ) break;

		save_GMM(savefile);
	}
	for(f=0;f<FileNum;f++){
		for(i=0;i<FeatLen[f];i++)
			free(likeli[f][i]);
		free(likeli[f]);
	}
	free(likeli);
	for(f=0;f<FileNum;f++)
		free(likeli_sum[f]);
	free(likeli_sum);
}



void K_means(void)
{  
	int ClusterNo,SplitWho,TotalSampleNum;
	int	i,f,t,j,k;
	int *SampleNum,**BelongTo;
	double	*dis;
	double TotalDistortion,e=0.02;

  	if((BelongTo=(int **)malloc(FileNum*sizeof(int *)))==NULL) Exit();
	for(f=0;f<FileNum;f++)
		if((BelongTo[f]=(int *)malloc(FeatLen[f]*sizeof(int)))==NULL) Exit();
	if((SampleNum=(int *)malloc(MixtureNum*sizeof(int)))==NULL) Exit();
	if((dis=(double *)malloc(MixtureNum*sizeof(double)))==NULL) Exit();

	for(f=0;f<FileNum;f++)
		for(t=0;t<FeatLen[f];t++)
			BelongTo[f][t]=0;

	/************************** non-binary spliting ***************************/
	ClusterNo=1;
	while(1){

    // calculate centroid
		for(k=0;k<ClusterNo;k++){
			for(j=0;j<FeatDim;j++)	GMMMean[k][j]=0.0;
			SampleNum[k]=0;
			for(f=0;f<FileNum;f++){
				for(t=0;t<FeatLen[f];t++){
					 if(BelongTo[f][t]==k){
						 for(j=0;j<FeatDim;j++)  GMMMean[k][j] +=(double)(Feat[f][t][j]);
						 SampleNum[k]++;
					 }
				}
			}
			if(SampleNum[k]==0) printf("%d null\n",k);
			for(j=0;j<FeatDim;j++) GMMMean[k][j] /=(double)SampleNum[k];
		}

		if(ClusterNo==MixtureNum) break;

		SplitWho=0;
		for(k=0;k<ClusterNo;k++)	if(SampleNum[k]>SampleNum[SplitWho]) SplitWho=k;
		for(j=0;j<FeatDim;j++){
			GMMMean[ClusterNo][j]=GMMMean[SplitWho][j]*(1.0+e);;
			GMMMean[SplitWho][j] *=(1.0-e);
		}
		ClusterNo++;

    // nearest neighbor search
		TotalDistortion=0.0;
		for(f=0;f<FileNum;f++){
			for(t=0;t<FeatLen[f];t++){
				for(k=0;k<ClusterNo;k++)  dis[k]=distance(Feat[f][t],GMMMean[k],FeatDim);
				BelongTo[f][t]=pick_minimum(dis,ClusterNo);
//      printf("%d %lf\n",BelongTo[i],dis[BelongTo[i]]);
				TotalDistortion +=dis[BelongTo[f][t]];
			}
		}
//		printf("%d %lf\n",ClusterNo,TotalDistortion);

	}
	/*********************************************************************/

    /** re-search and calculate centroids after cluster number has been fixed **/
	for(int re=0;re<5;re++){
	TotalDistortion=0.0;
	for(f=0;f<FileNum;f++){
		for(t=0;t<FeatLen[f];t++){
			for(k=0;k<ClusterNo;k++)  dis[k]=distance(Feat[f][t],GMMMean[k],FeatDim);
			BelongTo[f][t]=pick_minimum(dis,ClusterNo);
//      printf("%d %lf\n",BelongTo[i],dis[BelongTo[i]]);
			TotalDistortion +=dis[BelongTo[f][t]];
		}
	}
//	printf("%d %lf\n",ClusterNo,TotalDistortion);

	for(k=0;k<ClusterNo;k++){
		for(j=0;j<FeatDim;j++)	GMMMean[k][j]=0.0;
		SampleNum[k]=0;
		for(f=0;f<FileNum;f++){
			for(t=0;t<FeatLen[f];t++){
				 if(BelongTo[f][t]==k){
					 for(j=0;j<FeatDim;j++)  GMMMean[k][j] +=(double)(Feat[f][t][j]);
					 SampleNum[k]++;
				 }
			}
		}
		if(SampleNum[k]==0) printf("%d null\n",k);
		for(j=0;j<FeatDim;j++) GMMMean[k][j] /=(double)SampleNum[k];
	}
	/*************************************************************************/

  // calculate variance
	TotalSampleNum=0;
	for(k=0;k<MixtureNum;k++){
		for(j=0;j<FeatDim;j++)  GMMVar[k][j]=0.0;
		SampleNum[k]=0;
		for(f=0;f<FileNum;f++){
			for(t=0;t<FeatLen[f];t++){
				if(BelongTo[f][t]==k){
					for(j=0;j<FeatDim;j++)  GMMVar[k][j] +=((double)(Feat[f][t][j])-GMMMean[k][j])*((double)(Feat[f][t][j])-GMMMean[k][j]);
					SampleNum[k]++;
				}
			}
		}
		if(SampleNum[k]==0) printf("%d null\n",k);
		for(j=0;j<FeatDim;j++) GMMVar[k][j] /=(double)SampleNum[k];
		TotalSampleNum +=SampleNum[k];
	}

	for(i=0;i<MixtureNum;i++){
		for(j=0;j<FeatDim;j++){
			if(GMMVar[i][j]==0.0){
//				printf("k-means: covariance=0\n");
				GMMVar[i][j]=0.00000001;
			}
			GMMVarconst[i][j]=log(1.0/sqrt(GMMVar[i][j]*6.28318));
		}
	}

	for(k=0;k<MixtureNum;k++) GMMWeight[k]=(double)SampleNum[k]/(double)TotalSampleNum;
	}

  	free(dis);
	free(SampleNum);
	for(f=0;f<FileNum;f++) free(BelongTo[f]);
	free(BelongTo);
}


void save_GMM(char *safile)
{
	FILE *fp;
	int j;

	fp=fopen(safile,"wb");
	fwrite(&FeatDim,sizeof(int),1,fp);
	fwrite(&MixtureNum,sizeof(int),1,fp);
	for(j=0;j<MixtureNum;j++)
		fwrite(GMMMean[j],sizeof(double),FeatDim,fp);
	for(j=0;j<MixtureNum;j++)
		fwrite(GMMVar[j],sizeof(double),FeatDim,fp);
	fwrite(GMMWeight,sizeof(double),MixtureNum,fp);
//	printf("  saved!\n");
	fclose(fp);
}


void Exit(void)
{
  printf("memory insufficient\n");
  exit(0);
}


double GaussianPDF(float *GX,double *MeAn,double *VaR,double *VaRcOnSt,int Dim)
{
  double p=0.0;
  int i;

  for(i=0;i<Dim;i++)
    p += VaRcOnSt[i]+(((double)(GX[i])-MeAn[i])*((double)(GX[i])-MeAn[i])/(-2.0*VaR[i])) ;
  
  return p;
}

int pick_minimum(double *q,int numb)
{
  int i,j;
  double temp;

  temp=q[0];
  j=0;
  for(i=1;i<numb;i++){
	if(q[i]<=temp){
      j=i;
	  temp=q[i];
	}
  }

  return j;
}


double distance(float *A,double *B,int size)
{
  int i;
  double dis=0.0;
  for(i=0;i<size;i++)
	dis +=((double)(A[i])-B[i])*((double)(A[i])-B[i]);
  return dis;
}


void Free(void)
{
  int i,f;

  for(f=0;f<FileNum;f++){
    for(i=0;i<FeatLen[f];i++)
      free(Feat[f][i]);
    free(Feat[f]);
  }
  free(Feat);
  for(i=0;i<MixtureNum;i++){
      free(GMMMean[i]);
      free(GMMVar[i]);
      free(GMMVarconst[i]);
  }
  free(GMMMean);
  free(GMMVar);
  free(GMMVarconst);
  free(GMMWeight);

}
