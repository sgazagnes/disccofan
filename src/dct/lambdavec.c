#include "types.h"
#include "lambdavec.h"

/******************************************************************************/
/*                             Lambda Vector handling                         */
/******************************************************************************/

LambdaVec *lambda_vector_create(uint size)
{
   LambdaVec *lvec;

   lvec = malloc(sizeof(LambdaVec));
   if (lvec==NULL)  {
     warn("Allocation of lvec vector failed");
     return(NULL);
   }

   lvec->num_lambdas = size-1;
   lvec->lambdas = malloc(size*sizeof(float));
   if (lvec->lambdas==NULL)
   {
     warn("Allocation of lvec vector failed");
     free(lvec);
     return(NULL);
   }
   return(lvec);
}

void lambda_vector_delete(LambdaVec *lvec)
{
  free(lvec->lambdas);
  free(lvec);
} 

LambdaVec *lambda_vector_read(char *prefix, char *suffixe, double imScale)
{
   FILE *infile;
   LambdaVec *l;
   int size;
   int c;
   char* ts1 = strdup(prefix);
   char* dir = dirname(ts1);
   char *fname;
   asprintf(&fname,"%s/%s", dir, suffixe);
   infile = fopen(fname, "r");
   if (infile==NULL)  return(NULL);
   fscanf(infile, "LambdaVector\n");

   fscanf(infile, "%d\n", &size);

   l = lambda_vector_create(size);
   float test;
   if(l)
   {
     fscanf(infile,"%f ",&test);
     if(test != 0) error("Vector file must start with 0");
      for(c=0;c<size;c++)
	{   
	  fscanf(infile,"%f ",&l->lambdas[c]); 
	  l->lambdas[c] /= imScale;
	}  
   } 
   fclose(infile);
   return(l);
} 
