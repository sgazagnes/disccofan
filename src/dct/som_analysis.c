#include "types.h"
#include "som_analysis.h"
#include "flood.h"
/******************************************************************************/
/*                            SOM analysis                                    */
/******************************************************************************/

int *som_attribute_create(ulong size)
{
   int *som_attr;

   som_attr = (int *) calloc(size, sizeof(int));
   memset(som_attr, -1, size*sizeof(int));
   if (som_attr==NULL)  {
     warn("Allocation failed");
     return(NULL);
   }

   return(som_attr);
}


int *read_som_attributes(Arguments *args, ulong size,  char *prefix){
  FILE *infile;
  int size_som = args->somsize_arg;
  char* ts1 = strdup(prefix);
  char* dir = dirname(ts1);
  char *fname;
  asprintf(&fname,"%s/%s", dir, args->somfile_arg);
  infile = fopen(fname, "r");
  if (infile==NULL) {
    error("Couldn't read the lvec file at this address: %s", fname);
    return(NULL);
  }
  fscanf(infile, ",index,winning_Nx,winning_Ny\n");
  int *som_attr = som_attribute_create(size);
  int i, index, wNx, wNy;
  if(som_attr)
    {
      while(!feof(infile))
	{
	  fscanf(infile,"%d,%d,%d,%d",&i, &index, &wNx, &wNy);
     	  som_attr[index] = wNy*size_som + wNx;
	  //printf("%d \t ", som_attr[idx]);	 

	}  
    } 
  fclose(infile);
  return(som_attr);
}
      

void som_filter(Node *tree, value *out, int *som_attr, int lambda) {
  value  val;
  idx	 parent;
  ulong  v, u, w;
  bool *reached = calloc(tree->size_init, sizeof(bool));
  printf("HERE \n");
  for (v = 0; v < tree->size_init; v++) {//tree->size_init
    //   printf("%d, %d \n", v, reached[v]);
    if (!reached[v]) { /* not filtered yet */
      w =  v;
      parent = tree->parent[w];
      /* while ((parent != BOTTOM)) {
	ulong att = som_attr[w];
	if(som_attr[parent] != att){
	  printf("%d, %.9f, %.9f, %d, %d\n", w, tree->gval[w], tree->gval[parent],att, som_attr[parent]);
	}
        w = parent;
      parent = get_levelroot(tree,tree->parent[w]);
      }*/

      /* repeat while we're not at the bottom, th */
      while ((parent != BOTTOM)  && !reached[w] && ((tree->gval[w] == tree->gval[parent] )  ||  ( som_attr[w] != lambda))) {
        w = parent;
        parent = tree->parent[w];	      
      }
      // printf("OUT?\n");
      if (reached[w]) {
	val = out[w];
      } else if ( som_attr[w]== lambda) {
        val = tree->gval[w]; 
      } else {
        val = 0;
      }

      u = v;
      while (u != w) {
	//	if ((0 <= u) && (tree->size_init > u)){
	  out[u] 	= val;
	  reached[u] 	= true;
	  //	}
        u = tree->parent[u];
      }
      // if ( (lwb <= w) && (upb > w)){
	out[w] 		= val;
	reached[w]      = true;
	//	}
      }
  }
} /* tree_filter_direct */
