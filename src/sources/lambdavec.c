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

   lvec->num_lambdas = size;
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

LambdaVec *lambda_vector_read(char *name, double im_scale)
{
    FILE *infile;
    LambdaVec *l;
    int size = 0;
    int c;
    float value;
    
    infile = fopen(name, "r");
    if (infile == NULL) {
        error("Couldn't read the lvec file: %s. Make sure the path is correct", name);
        return NULL;
    }

    // Skip the first line
    fscanf(infile, "#LambdaVector file. Keep this exact line in all similar lambda vector files\n");

    // Count the number of lambda values
    long initial_position = ftell(infile); // Save the current file position
    while (fscanf(infile, "%f", &value) == 1) {
        size++;
    }

    // Create the lambda vector
    l = lambda_vector_create(size+1); //(We are adding a 0 value)
    if (l == NULL) {
        fclose(infile);
        return NULL;
    }

    // Rewind file to start reading lambdas
    fseek(infile, initial_position, SEEK_SET);

    // Read values and store them in lambda vector
    // fscanf(infile, "%f", &value);
    // if (value != 0) {
    //     error("Vector file must start with 0");
    //     fclose(infile);
    //     return NULL;
    // }

    l->lambdas[0] = 0;
    for (c = 0; c < size; c++) {
        fscanf(infile, "%f", &l->lambdas[c+1]);
        l->lambdas[c+1] /= im_scale;
    }

    fclose(infile);
    return l;
}

void print_lambda_vec(LambdaVec *lvec)
{
  info("Printing Lambda vector values read");
  info("Number of lambdas in vector: %d", lvec->num_lambdas-1);
  info("Lambdas:");
  if(rank() == 0){
    for (int i = 1; i < lvec->num_lambdas; i++)
      printf("%f\n", lvec->lambdas[i]);
  }
}
