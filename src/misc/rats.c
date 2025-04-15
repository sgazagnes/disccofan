
/****************************************************************************/
/*                                                                          */
/* module rats.c: implements RATS algorithm for local, bilevel              */ 
/*                thresholding, using square-sobel filtering                */
/*                                                                          */
/* References:    M.H.F. Wilkinson (1998) Optimizing edge detectors for     */
/*                robust automatic threshold selection. Graph. Models       */
/*                Image Proc. 60:                                           */
/*                M.H.F. Wilkinson (1996) Rapid automatic segmentation of   */
/*                fluorescent and phase-contrast images of bacteria. In:    */
/*                Fluorescence Microscopy and Fluorescent Probes,           */
/*                (J. Slavik, ed), pp 261-266, Plenum Press, New York.      */ 
/* Author:        Michael H. F. Wilkinson                                   */
/*                Institute for Mathematics and Computing Science           */
/*                University of Groningen,                                  */
/*                PO Box 800, 9700 AV Groningen, The Netherlands            */
/*                e-mail: michael@cs.rug.nl                                 */
/* Version:       14-04-2000                                                */
/* Comments:      Feel free to use non-commercially, but do acknowledge     */
/*                copyright, or (even better) cite articles above  :-)      */
/*                Contains main program for testing purposes (reads .pgm    */
/*                files, produces file rats.pgm                             */
/****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <FreeImage.h>

/**********************************************************************/
/*                                                                    */
/*                           Quadtree type                            */
/*                                                                    */
/**********************************************************************/


typedef double **StatImage;   
         /* 2D array to store each level of a quadtree */
 
typedef struct 
          { int numlevels,    /* number of levels in the quadtree        */ 
	        maxsubdivs,   /* maximum subdivisions (=2^(numlevels-1)) */
                xleafsize,    /* x-size of leaves                        */
                yleafsize;    /* y-size of leaves                        */
            StatImage *quadtree;  /* the quadtree itself                 */
          } StatPyramid;

/* NOTE that zero is the root level, numlevels-1 the leaf level */
/* slightly confusingly, I refer to the 0 level as the highest  */
/* which makes sense if you think of pyramids                   */

/**********************************************************************/
/*                                                                    */
/*                   Handlers for quadtree type                       */
/*                                                                    */
/**********************************************************************/

void InitStatPyramid ( StatPyramid *pyramid,
                       int levels,
                       int imagewidth,
                       int imageheight      )
/* Initializes StatPyramid type creating the quadtree and computing the
   leaf sizes 
*/
{ unsigned int 
    i,j,k,      /* indices               */
    max=1;      /* size of current level */

  /* initialize number of levels */
  pyramid->numlevels=levels;      

  /* allocate memory for pointers to pyramid levels */   
  pyramid->quadtree=(StatImage *)malloc(levels*sizeof(StatImage));
 
  /* for all levels do */
  for (i=0;i<levels;i++,max*=2)
    { /* allocate column of pointers to StatImage rows */ 
      pyramid->quadtree[i]= (StatImage)malloc(max * sizeof(double *));
      /* for all rows pointers */
      for(j=0;j<max;j++)
        { /* allocate rows */
	  pyramid->quadtree[i][j]=(double *)malloc(max*sizeof(double));
          /* for all elements of row set to zero */
          for(k=0;k<max;k++)
            pyramid->quadtree[i][j][k]=0.0;
	}
    }
  /* compute maximum number of subdivisions at leaf level */
  max /=2;
  pyramid->maxsubdivs=max;

  /* compute x and y leaf dimensions */
  pyramid->xleafsize = (imagewidth+max -1) / max;
  pyramid->yleafsize = (imageheight+max-1) / max;

}

void ExitStatPyramid ( StatPyramid *pyramid )
/* Frees contents of StatPyramid type, and clears levels, etc */
{ int i,j,max=1;
  for (i=0;i<pyramid->numlevels;i++,max*=2)
    { for (j=0;j<max;j++)
	free(pyramid->quadtree[i][j]);       
      free(pyramid->quadtree[i]);
    }
  free(pyramid->quadtree);
  pyramid->numlevels=0;   
}

/**********************************************************************/
/*                                                                    */
/* Pixel type: can be changed to other type for which comparison with */
/* double type is defined                                             */
/*                                                                    */
/**********************************************************************/

typedef double Pixel;

/**********************************************************************/
/*                                                                    */
/*                           2D image type                            */
/*                                                                    */
/**********************************************************************/

typedef struct 
          { int width,     /* image width  */
                height;    /* image height */
            Pixel **data;  /* pointer to 2D array */
          } Image2D;

/**********************************************************************/
/*                                                                    */
/*                    Handlers for Image2D type                       */
/*                                                                    */
/**********************************************************************/

void InitImage2D ( Image2D *image,
                   int     width,
                   int     height )
/* Initialize image to width * height  */
{ int y;
  image->width=width;
  image->height=height;
  image->data= (Pixel **)malloc(height * sizeof(Pixel *));
  for (y=0;y<height;y++)
     image->data[y]=(Pixel *)malloc(width*sizeof(Pixel));
}

void ExitImage2D ( Image2D *image )
/* Clears contents of image, freeing memory */
{ int y;
  for (y=0;y<image->height;y++)
    free(image->data[y]);
  free(image->data);
  image->height=0;
  image->width=0;
}

void RowRead(int x, int y, Image2D *image, int width, Pixel *buf)
/* Reads part of image into the pixel buffer buf, */
/* starting at x,y, of length width               */
{ int i;
  for (i=0;i<width;i++)
    buf[i]=image->data[y][x+i]; 
}

/**********************************************************************/
/*                                                                    */
/*                    Quadtree filling routines                       */
/*                                                                    */
/**********************************************************************/
 
void FillLowest( StatPyramid *sums, 
                 StatPyramid *weights,
                 Image2D     *in,
                 double      lambdasigma   )
/* Fills lowest level of both quadtrees denominators (weights) and */
/* enumerators (sums) of RATS' T-statistic                         */
{ int x,y;    /* take a wild guess */
  int k;
  double 
    *cursum,    /* pointer to current leaf in sums quadtree*/ 
    *curweight, /* pointer to current leaf in weights quadtree*/ 
    weight,     /* pixel weight */
    wx,         /* weight contribution of x-gradient */ 
    wy;         /* weight contribution of y-gradient */ 
  /* Pixel buffers */
  Pixel *temp,                                 /* temporary for swapping */
    *prev=malloc((in->width+2)*sizeof(Pixel)), /* previous Pixel row     */
    *cur= malloc((in->width+2)*sizeof(Pixel)), /* current Pixel row      */
    *next=malloc((in->width+2)*sizeof(Pixel)); /* next Pixel row         */

  /* increment pixel buffers so the array is indexed from -1 to width    */
  prev++; 
  next++;
  cur++;

  /* read initial buffers */
  RowRead(0,0,in,in->width,prev);
  RowRead(0,0,in,in->width,cur);

  /* copy image edge values into Pixels adjacent to image */
  prev[-1]=prev[0];
  prev[in->width]=prev[in->width-1];

  cur[-1]=cur[0];
  cur[in->width]=cur[in->width-1];

  /* scan entire image */
  for (y=0;y<in->height;y++)
    { /* read next row (y+1), unless y=height-1, then read row y */
      RowRead(0,y+(y!=in->height-1),in,in->width,next);

      /* init leaf pointers */
      cursum = sums->quadtree[weights->numlevels-1][y/sums->yleafsize];
      curweight = weights->quadtree[weights->numlevels-1][y/sums->yleafsize];

      /* pad buffer as above */
      next[-1]=next[0];
      next[in->width]=next[in->width-1];

      /* scan image row */
      
      for (x=0, k=0; x<in->width;x++, k++)
        { 
          if (k==sums->xleafsize)
	    /* if boundary of leaf is reached */
	    { /* move to next leaves */
              cursum ++;
              curweight ++;
              /* reset within leaf counter */
	      k=0;
	    }
          /* Sobel gradient in x direction */ 
          wx=(double)(prev[x-1])-(double)(next[x-1])+
             2*((double)(prev[x])-(double)(next[x]))+
             (double)(prev[x+1])-(double)(next[x+1]);

	  /* Sobel Gradient in y direction */
          wy=(double)(prev[x-1])-(double)(prev[x+1])+
             2*((double)(cur[x-1])-(double)(cur[x+1]))+
             (double)(next[x-1])-(double)(next[x+1]);
          /* compute square Sobel gradient (normalized) */
          weight=(wx*wx + wy*wy)/16.0;
          /* check whether weight is large enough */
          if (weight>lambdasigma)
            { /* add weight to correct leaf of weights */
              (*curweight) +=weight;
              /* add weight * current greylevel to correct leaf of sums */
              (*cursum) +=weight*(double)(cur[x]);
	    }
        }
      /* swap Pixel row buffers */
      temp=prev;prev=cur;cur=next;next=temp;
    }
  /* clean up Pixel row buffers (not forgetting to decrement them first) */ 
  next--;
  free(next);
  prev--;
  free(prev);
  cur--;
  free(cur);
}

void FillOthers (StatPyramid *stat)
/* Fill all other levels of the quadtree in stat based on leaf values */
{ unsigned int 
    i,j,k,                /* quadtree indices */
    max=stat->maxsubdivs; /* size of leaf level */

  /* start at leaf level, moving to root */
  for (i=(stat->numlevels-1);i>0;i--,max/=2)
    /* for all nodes at each level */
    for(j=0;j<max;j++)
      for(k=0;k<max;k++)
        /* add current nodes at current level to those one higher */
        (stat->quadtree[i-1][j/2][k/2])+= (stat->quadtree[i][j][k]);
}


/**********************************************************************/
/*                                                                    */
/*           type containing thresholds of leaf centroids             */
/*                                                                    */
/**********************************************************************/

typedef double **ThreshMap;


/**********************************************************************/
/*                                                                    */
/*                   Handlers for ThreshMap type                      */
/*                                                                    */
/**********************************************************************/


ThreshMap CreateThreshMap(int size )
/* allocates memory for size*size array of double */
{ int i;
  ThreshMap newmap= (ThreshMap) malloc(size*sizeof(double *));
  for (i=0;i<size;i++)
    newmap[i]=(double *) malloc(size*sizeof(double));
  return newmap;
}

void FreeThreshMap (ThreshMap thresh, int size )
/* frees the contents of thresh              */
/* size must be identical to that allocated!!*/
{ int i;
  for (i=0;i<size;i++)
    free(thresh[i]);
  free(thresh);
}

/**********************************************************************/
/*                                                                    */
/*            Thresholding routines proper (at last)                  */
/*                                                                    */
/**********************************************************************/

double RatsThresh (StatPyramid *sums,  
                   StatPyramid *weights,
                   int i, int j, int k,
                   double lambdasigma   )
/* Recursive routine computing rats threshold for centroids of leaf */
/* indicated by j and k. To be called with i=sums->numlevels, and   */
/* lambdasigma as above.                                            */
{ double retval;
  /* if current denominator indicates presence of edge */
  if ((weights->quadtree[i][j][k]>3*lambdasigma) || (i==0) )
    /* use T-statistic at current level */
    retval=(weights->quadtree[i][j][k]>0) ? 
           sums->quadtree[i][j][k]/weights->quadtree[i][j][k]: 128;
  else 
    /* move one up in the hierarchy */
    retval=RatsThresh(sums,weights,i-1,j/2,k/2,lambdasigma);
  return retval;
}


void ThreshAll ( Image2D *in, Image2D *out, 
                 ThreshMap thresh, StatPyramid *stat )
/* Performs bilinear interpolation between the leaf-centroids   */
/* surrounding all (x,y) to compute a local threshold, computes */
/* thresholded image                                            */
{ int 
    x,y,        /* x and y indices */
    j0,         /* index of lefthand leaves  */
    j1,         /* index of righthand leaves */
    k0=0,       /* index of top leaves       */
    k1=0,       /* index of bottom leaves    */
    countx,     /* within leaf pixel x index */
    county;     /* within leaf pixel y index */

  double 
    val0, val1,    /* interim threshold values after x interpolation */
    wx1,           /* weights of righthand leaf thresholds           */
    wy1;           /* weights of bottom leaf thresholds              */

  county=stat->yleafsize/2;             /* initialize halfway through leaf */
  for (y=0;y<in->height;y++, county++)
    { if (county==stat->yleafsize)                       /* if end of leaf */
        { k0+=(k1>0);                    /* increment k0 UNLESS k1 is zero */
	  k1+=(k0<stat->maxsubdivs-1);   /* increment k1 UNLESS k0 is max  */
          county=0;                      /* reset counter                  */
	}

      /* compute y-weighing factor */
      wy1=((double)county-0.5) /(double)stat->yleafsize;      
      
      /* initialy, left and right leaf indexes are zero */
      j1=0; 
      j0=0;
      countx=stat->xleafsize/2;      /* initialize halfway through leaf */
      for (x=0;x<in->width;x++, countx++)
        { if (countx==stat->xleafsize)                    /* if end of leaf */
            { j0+=(j1>0);                 /* increment j0 UNLESS j1 is zero */
              j1+=(j0<stat->maxsubdivs-1); /* increment j1 UNLESS j0 is max */
              countx=0;                    /* reset counter                 */
	    }
          
          /* compute x-weighing factor */
          wx1=((double)countx-0.5) /(double)stat->xleafsize;

          /* perform x-interpolation */
          val0=wx1*thresh[k0][j1]+(1-wx1)*thresh[k0][j0];
          val1=wx1*thresh[k1][j1]+(1-wx1)*thresh[k1][j0];

          /* perform y-interpolation and threshold */
          out->data[y][x]=in->data[y][x]>(wy1*val1+(1-wy1)*val0) ? 255 : 0;
	  
     	}
    } 
}

void RatsImage2D ( Image2D *in, Image2D *out, 
                   int levels, double  noise, double  lambda )  
/* RatsImage2D yields a thresholded image, allowing free setting */
/* of number of levels in quadtree, and noise and lambda levels  */
/* parameters:                                                   */
/*         in      :      input image                            */
/*         out     :      output image (must be same size as in) */
/*         levels  :      levels in quadtree                     */
/*         noise   :      noise level (standard deviation)       */
/*         lambda  :      pixels with gradients below            */
/*                        lambda*noise are not used in           */
/*                        computation of threshold.              */
/*                        recommended value: 3.0                 */ 
{ unsigned int j,k,x,y;             /* various indices */
  StatPyramid sums,weights;         /* quadtrees       */
  double lambdasigma=lambda*noise;  /* edge threshold value */
  ThreshMap thresh;                 /* thresholds of leaves */

  /* initialize auxiliary data types */
  InitStatPyramid(&sums,levels,in->width,in->height);
  InitStatPyramid(&weights,levels,in->width,in->height);

  thresh=CreateThreshMap(sums.maxsubdivs);

  /* square gradient threshold to account for use of square of gradients */
  lambdasigma*=lambdasigma;         


  /* fill quadtree levels */
  FillLowest (&sums,&weights,in,lambdasigma );

  FillOthers(&sums);
  FillOthers(&weights);



  /* compute thresholds for leaf centroids */
  for (j=0;j<sums.maxsubdivs;j++)
    for (k=0;k<sums.maxsubdivs;k++)
      thresh[j][k]=RatsThresh(&sums,&weights,levels-1,j,k,lambdasigma);

  /* the actual thresholding */
  ThreshAll(in,out,thresh,&sums);
/*  for (y=0;y<in->height;y++)
      for (x=0;x<in->width;x++)
        out->data[y][x]=(in->data[y][x]>ThreshVal(&sums,thresh,x,y)) ? 255 :0;
*/
  /* clean up auxiliary data */
  FreeThreshMap(thresh,sums.maxsubdivs);
  
  ExitStatPyramid(&sums);
  ExitStatPyramid(&weights);

}


/***************************************************************************/
/*                                                                         */
/*          Some quick and dirty PGM-file import and export                */
/*                                                                         */
/***************************************************************************/
FIBITMAP* GenericLoader(const char* lpszPathName, int flag) {
  FREE_IMAGE_FORMAT fif = FIF_UNKNOWN;
//  check the file signature and deduce its format
// (the second argument is currently not used by FreeImage)
  fif = FreeImage_GetFileType(lpszPathName, 0);
  if(fif == FIF_UNKNOWN) {
    // no signature ?
    // try to guess the file format from the file extension
    fif = FreeImage_GetFIFFromFilename(lpszPathName);
  }
  // check that the plugin has reading capabilities ...
  if((fif != FIF_UNKNOWN) && FreeImage_FIFSupportsReading(fif)) {
    // ok, let's load the file
    FIBITMAP *dib = FreeImage_Load(fif, lpszPathName, flag);
    // unless a bad file format, we are done !
    return dib;
  }
  return NULL;
}

unsigned long ReadTIFF(char *fnm, Image2D *im){
   FIBITMAP *dib = GenericLoader(fnm,0);
   unsigned long  bitsperpixel;
   unsigned int x,y,i,imsize,numplanes=1;
   if (dib == NULL) return 0;
     
   bitsperpixel =  FreeImage_GetBPP(dib);
   if ((bitsperpixel==24) ||(bitsperpixel==48)){
     numplanes=3;
     bitsperpixel/=3;
   }
   InitImage2D(im, FreeImage_GetWidth(dib),FreeImage_GetHeight(dib));

   printf("BitsPerPixel = %d, Width= %d, Height= %d\n", bitsperpixel, im->width, im->height);
   switch(bitsperpixel) {
   case 8:
     
     for(y = 0; y < im->height; y++) {
       BYTE *bits = (BYTE *)FreeImage_GetScanLine(dib, y);
       for(x = 0; x < im->width; x++,i++) {
	 im->data[y][x] = bits[numplanes*x];
       }
     }
     
     FreeImage_Unload(dib);
     break;
   case 16:
     i=0;
     for(y = 0; y < im->height; y++) {
       unsigned short *bits = (unsigned short *)FreeImage_GetScanLine(dib, y);
       for(x = 0; x < im->width; x++,i++) {
	 im->data[y][x] = bits[numplanes*x];
       }
     }
     FreeImage_Unload(dib);
     break; 
   default : 
     FreeImage_Unload(dib);
     
      
   }
   return bitsperpixel;  
}

void WriteTIFF( char *fname, Image2D *im, unsigned long bitspp){
  FIBITMAP *outmap;
  long i,j,y,x; 
  FREE_IMAGE_FORMAT fif = FreeImage_GetFIFFromFilename(fname);
  int maxgrey;
  if (bitspp == 8){
    unsigned char *imagebuf;
    RGBQUAD *pal;
    maxgrey=255;
    outmap = FreeImage_AllocateT(FIT_BITMAP,im->width,im->height,bitspp,0xFF,0xFF,0xFF);
    pal = FreeImage_GetPalette(outmap);
    for (i = 0; i < 256; i++) {
      pal[i].rgbRed = i;
      pal[i].rgbGreen = i;
      pal[i].rgbBlue = i;
    }
    i = 0;
    for (y=0; y< im->height; y++){      
      imagebuf = FreeImage_GetScanLine(outmap,y);
      for (x=0;x<im->width;x++)
	imagebuf[x]=(im->data[y][x]<=maxgrey)?im->data[y][x]:maxgrey;
	
    }

  } else {
    unsigned short *imagebuf;
    outmap = FreeImage_AllocateT(FIT_UINT16,im->width,im->height,16,0xFFFF,0xFFFF,0xFFFF);
     maxgrey=0xFFFF;
   for (y=0; y<im->height; y++){      
      imagebuf = (unsigned short *)FreeImage_GetScanLine(outmap,y);
      for (x=0;x<im->width;x++)
	 imagebuf[x]=(im->data[y][x]<=maxgrey)?im->data[y][x]:maxgrey;
	
    }
  }
  FreeImage_Save(fif,outmap,fname,0); 
  FreeImage_Unload(outmap);

}




/***************************************************************************/
/*                                                                         */
/*                   Main routine for testing purposes                     */
/*                                                                         */
/***************************************************************************/


int main (int argc, char **argv)
{
  int levels=5;
  double noise = 2;
  Image2D in, out;
  unsigned long bpp;

  bpp = ReadTIFF(argv[1], &in );
  
  InitImage2D(&out, in.width, in.height);

  if (argc>2) levels=atoi(argv[2]);
  if (argc>3) noise=atof(argv[3]);

  #RatsImage2D( &in, &out, levels, noise, 3.0 );

  WriteTIFF("rats.tif", &in, 8 );
  ExitImage2D(&in);
  ExitImage2D(&out);

  return;
}







