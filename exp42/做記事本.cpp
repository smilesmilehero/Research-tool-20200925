#include "stdio.h"
#include <stdlib.h>
int main(int argc, const char * argv[])
{

FILE *fp;
 
int k;


fp = fopen("k.txt", "w");//檔名 
if(fp ==NULL)
printf("打開文件失敗\n");

		for(k=1;k<11;k++)
 	

fprintf(fp,"gmmdrec tedt.crp ref/GMM-%d.16 16 %d   40 \n",k,k);					           

 return 0;

} 



