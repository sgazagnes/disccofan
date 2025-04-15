/**************************************************************************/
/*                                                                        */
/*                                                                        */
/*      recgauss.c : recursive Gaussian filter using method of Young      */
/*                   and Van Vliet (1995 Signal Proc. 44:139-151)         */
/*                                                                        */
/*                                                                        */
/*                                                                        */
/*                                                                        */
/*                                                                        */
/*                                                                        */
/*                                                                        */
/*                                                                        */
/*                                                                        */
/**************************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <malloc.h>
#include <time.h>
#include <sys/times.h>
#include <unistd.h>
#include <FreeImage.h>

#define PI 3.1415926


void fillGaussKernel(double sigma, double *kernel, int size)
{ 
  int i;
  for (i=0;i<size;i++)
    kernel[i]=exp(-(double)(i*i) / sigma)/(sigma*sqrt(2*PI));
}

void bluntGauss2DFilter( double **in, 
			 double **out, 
			 int width, 
			 int height, 
			 double sigma )
{ 
  int x,y,i, size=8*(int)(sigma+0.5);
  double *kernel = (double *)malloc(size * sizeof(double));
  double *buf = (double *)malloc(height * sizeof(double));
  double sum,weight;
  fillGaussKernel(sigma,kernel,size);
  for(y=0;y<height;y++)
    for(x=0;x<width;x++)
      {
	sum=0.0;
	weight=0.0;
	for (i=-size+1;i<size;i++)
	  if (((x+i)>=0) && ((x+i)<width))
	    { 
	      sum += kernel[abs(i)]*in[y][x+i];
	      weight += kernel[abs(i)];
	    }
	out[y][x]=sum/weight;
      }
  for(x=0;x<width;x++)
    {
      for(y=0;y<height;y++)
	{
	  sum=0.0;
	  weight=0.0;
	  for (i=-size+1;i<size;i++)
	    if (((y+i)>=0) && ((y+i)<height))
	      { 
		sum += kernel[abs(i)]*out[y+i][x];
		weight += kernel[abs(i)];
	      }
	  buf[y]=sum/weight;
	  
	}
      for(y=0;y<height;y++)
	out[y][x]=buf[y];
    }
  free(kernel);
  free(buf);
}

void computeBvalues( double *B, double *b0, double *b1, 
		     double *b2, double *b3, double q )
{
  *b0 = 1.57825 + 2.44413 * q + 1.4281 * q * q + 0.422205 * q * q * q;

  *b1 =( 2.44413 * q + 2.85619 * q * q + 1.26661  * q * q * q ) / *b0;

  *b2 =  -( 1.4281 * q * q + 1.26661  * q * q * q ) / *b0;

  *b3 =  (0.422205 * q * q * q)/ *b0;
   
  *b0 = 1;

  *B = 1.0  - ( *b1 + *b2 + *b3 );
}


double computeQvalue( double sigma )
{ 
  if (sigma > 2.5)
    return ( 0.98711*sigma - 0.96330 );
  else if (sigma > 0.5)
    return ( 3.97156 - 4.14554*sqrt(1-0.26891*sigma) );
  else
    return 0;   /* method not defined for sigma <0.5 */
}

void forwardFilter ( double *in, double *w, int len, double B,
                     double b0, double b1, double b2, double b3 )
{  
  int i;
  w[0] = in[0];
  w[1] = B*in[1] + (1 - B)*w[0];
  w[2] = B*in[2] + ( b1 * w[1] + (b2+b3) * w[0] );
  
  for (i = 3; i<len; i++)
    {
      w[i] = B * in[i] + b1 * w[i-1] + b2 * w[i-2] +b3 * w[i-3];
    }
}  

void backwardFilter ( double *w, double *out, int len, double B,
                     double b0, double b1, double b2, double b3 )
{  
  int i;
  out[len-1] = w[len-1];
  out[len-2] = B*w[len-2] + (1 - B)*out[len-1];
  out[len-3] = B*w[len-3] + ( b1 * out[len-2] + (b2+b3) * out[len-1] );
  
  for (i = len-4; i >= 0 ; i--)
    {
      out[i] = B * w[i] + ( b1 * out[i+1] + b2 * out[i+2] +b3 * out[i+3] );
    }
}  

void forward2DYFilter ( double **in, double **w, 
                        int width, int height, 
                        double B, double b0, 
                        double b1, double b2, double b3 )
{  
  int x, y;
  for (x=0;x<width;x++)
    { 
      w[0][x] = in[0][x];
      w[1][x] = B*in[1][x] + (1 - B)*w[0][x];
      w[2][x] = B*in[2][x] + ( b1 * w[1][x] + (b2+b3) * w[0][x] );
    }
   
   for (y = 3; y<height; y++)
     for (x=0;x<width;x++)
       {
	  w[y][x] = B * in[y][x] + 
	    ( b1 * w[y-1][x] + b2 * w[y-2][x] + b3 * w[y-3][x] );
       }
}  

void backward2DYFilter ( double **w, double **out, int width, int height,  
                         double B,
                         double b0, double b1, double b2, double b3 )
{  
  int x,y;
  for (x=0;x<width;x++)
    { 
      out[height-1][x] = w[height-1][x];
      out[height-2][x] = B*w[height-2][x] + (1 - B)*out[height-1][x];
      out[height-3][x] = B*w[height-3][x] + 
	( b1 * out[height-2][x] + (b2+b3) * out[height-1][x] );
    }
  for (y = height-4; y >= 0 ; y--)
     for (x=0;x<width;x++)
       {
	 out[y][x] = B * w[y][x] + 
	   ( b1 * out[y+1][x] + b2 * out[y+2][x] +b3 * out[y+3][x] );
       }
} 
 
void recGauss2DFilter( double **in, 
		       double **out, 
                       int width, 
		       int height, 
		       double sigma )
{ 
  double B, b0, b1, b2, b3, q = computeQvalue(sigma); 
  int x,y;
  computeBvalues(&B, &b0, &b1, &b2, &b3, q);
  
  for(y=0;y<height;y++)
    {
      forwardFilter(in[y], out[y], width, B, b0, b1, b2, b3);
      backwardFilter(out[y], out[y], width, B, b0, b1, b2, b3);
    }
  forward2DYFilter(out,out,width,height,B, b0, b1, b2, b3);
  backward2DYFilter(out,out,width,height,B, b0, b1, b2, b3);
}


void recGaussFilter( double *in, 
		     double *out, 
                     int len, 
		     double sigma )
{ 

  if (sigma<=50.0 && ((sigma/(double)len)<0.1 ) )
    { 
      double B, b0, b1, b2, b3, q = computeQvalue(sigma); 
      computeBvalues(&B, &b0, &b1, &b2, &b3, q);
      
      forwardFilter(in, out, len, B, b0, b1, b2, b3);
      backwardFilter(out, out, len, B, b0, b1, b2, b3);
    }
  else
    {
      recGaussFilter(in,out,len,sigma/2.0);
      recGaussFilter(out,out,len,sigma/2.0);
      recGaussFilter(out,out,len,sigma/2.0);
      recGaussFilter(out,out,len,sigma/2.0);
    }
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


/***************************************************************************/
/*                                                                         */
/*                          Some Point Operators                           */
/*                                                                         */
/***************************************************************************/


void mulImage2D(Image2D *in1, Image2D *in2, Image2D *out)
{
  int x,y;

  for (y=0;y<in1->height;y++)
    {
      for (x=0;x<in1->width;x++)
	out->data[y][x]=in1->data[y][x]*in2->data[y][x];
    }
}

void setImage2D(Image2D *out, Pixel val)
{
  int x,y;

  for (y=0;y<out->height;y++)
    {
      for (x=0;x<out->width;x++)
	out->data[y][x]=val;
    }
}

Pixel sumGreyImage2D(Image2D *out)
{
  int x,y;
  Pixel val=0.0;

  for (y=0;y<out->height;y++)
    {
      for (x=0;x<out->width;x++)
	val += out->data[y][x];
    }
  return val;
}

void subtractImage2D(Image2D *in1, Image2D *in2, Image2D *out)
{
  int x,y;

  for (y=0;y<in1->height;y++)
    {
      for (x=0;x<in1->width;x++)
	out->data[y][x]=(in1->data[y][x]-in2->data[y][x]+255)/2.0 ;
    }
}



void divImage2D(Image2D *in1, Image2D *in2, Image2D *out, 
		Pixel lowlimit, Pixel defthresh)
{
  int x,y;

  for (y=0;y<in1->height;y++)
    {
      for (x=0;x<in1->width;x++)
	out->data[y][x]=(in2->data[y][x]>lowlimit) ? 
	  in1->data[y][x]/in2->data[y][x] : defthresh; 
    }
}

void clipLowImage2D(Image2D *in, Image2D *out, Pixel lambdasigma)
{
  int x,y;

  for (y=0;y<in->height;y++)
    {
      for (x=0;x<in->width;x++)
	out->data[y][x]=(in->data[y][x]>lambdasigma)?in->data[y][x]:0;
    }
}

void threshImage2D(Image2D *in1, Image2D *in2, Image2D *out)
{
  int x,y;

  for (y=0;y<in1->height;y++)
    {
      for (x=0;x<in1->width;x++)
	out->data[y][x]=255*(in1->data[y][x]>in2->data[y][x]);
    }
}



/***************************************************************************/
/*                                                                         */
/*                         Squared Gradient filter                         */
/*                                                                         */
/***************************************************************************/
#define SQR(x) ((x)*(x))

void RowRead(int x, int y, Image2D *image, int width, Pixel *buf)
/* Reads part of image into the pixel buffer buf, */
/* starting at x,y, of length width               */
{ int i;
  for (i=0;i<width;i++)
    buf[i]=image->data[y][x+i]; 
}

void SquareSobel( Image2D     *in, Image2D *out )                 
{ int x,y;    /* take a wild guess */
  int k;
  double 
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

      /* pad buffer as above */
      next[-1]=next[0];
      next[in->width]=next[in->width-1];

      /* scan image row */
      
      for (x=0, k=0; x<in->width;x++, k++)
        { 
          /* Sobel gradient in x direction */ 
          wx=(double)(prev[x-1])-(double)(next[x-1])+
             2*((double)(prev[x])-(double)(next[x]))+
             (double)(prev[x+1])-(double)(next[x+1]);

	  /* Sobel Gradient in y direction */
          wy=(double)(prev[x-1])-(double)(prev[x+1])+
             2*((double)(cur[x-1])-(double)(cur[x+1]))+
             (double)(next[x-1])-(double)(next[x+1]);
          /* compute square Sobel gradient (normalized) */
          out->data[y][x]=(wx*wx + wy*wy)/16.0;
	  

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

void squaredGradient(Image2D *in, Image2D *out)
{
  int x,y;
  /* first row */  
  out->data[0][0]=SQR(in->data[0][0] - in->data[0][1]) +
    SQR(in->data[0][0] - in->data[1][0]);
  for(x=1;x<in->width-1;x++)
    {
      out->data[0][x]=SQR(in->data[0][x-1] - in->data[0][x+1]) +
        SQR(in->data[0][x] - in->data[1][x]);
      
    }
  out->data[0][in->width-1]=
    SQR(in->data[0][in->width-2] - in->data[0][in->width-1]) +
    SQR(in->data[0][in->width-1] - in->data[1][in->width-1]);
  
  /* all but first and last rows */

  for (y=1; y<in->height-1;y++)
    {
      out->data[y][0]=SQR(in->data[y][0] - in->data[y][1]) +
	SQR(in->data[y-1][0] - in->data[y+1][0]);
      for(x=1;x<in->width-1;x++)
	{
	  out->data[y][x]=SQR(in->data[y][x-1] - in->data[y][x+1]) +
	    SQR(in->data[y-1][x] - in->data[y+1][x]);
	  
	}
      out->data[0][in->width-1]=
	SQR(in->data[y][in->width-2] - in->data[y][in->width-1]) +
	SQR(in->data[y-1][in->width-1] - in->data[y+1][in->width-1]);
      
    }

  /* last row */
  out->data[in->height-1][0]=
    SQR(in->data[in->height-1][0] - in->data[in->height-1][1]) +
    SQR(in->data[in->height-2][0] - in->data[in->height-1][0]);
  for(x=1;x<in->width-1;x++)
    {
      out->data[in->height-1][x]=
	SQR(in->data[in->height-1][x-1] - in->data[in->height-1][x+1]) +
        SQR(in ->data[in->height-2][x] - in->data[in->height-1][x]);
      
    }
  out->data[in->height-1][in->width-1]=
    SQR(in->data[in->height-1][in->width-2] 
	- in->data[in->height-1][in->width-1]) +
    SQR(in ->data[in->height-2][in->width-1] 
	- in->data[in->height-1][in->width-1]);
}

void HistogramImage2D (Image2D *in, int maxgrey, long *hist)
{ 
  int x,y,j;
  for (j=0;j<maxgrey;j++)
    hist[j]=0;

  for (y=0; y<in->height;y++)
    {
      for(x=0;x<in->width;x++)
	{
          if (in->data[y][x]<maxgrey)
	    hist[(int)in->data[y][x]]++;
	  
	}
      
    }

}


float gasdev(int *idum);


/***************************************************************************/
/*                                                                         */
/*                         Illumination Functions                          */
/*                                                                         */
/***************************************************************************/

void slopeIlluminateImage2D(Image2D *in, Image2D *out,
                            float leftlevel, float rightlevel)
{ int x,y;

 for (y=0;y<in->height;y++)
   {
     for (x=0;x<in->width;x++)
       out->data[y][x]= in->data[y][x] * 
	 (rightlevel*x +leftlevel*(in->width-1-x))/((float)in->width-1);
   }
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



double  T_Gw(double lambda, double sigma, double eta)
{
  return 4* eta * eta*exp(-lambda*lambda/8.0) * (1 + lambda*lambda/4.0)
    * (exp(-lambda*lambda/8.0)+3/(2*sigma*sqrt(PI)));
}



#define BIG 1e36



void GM_RATS(Image2D *in, Image2D *out, 
	     double lambda, double eta, 
	     double sigmin, double sigmax )
{
  Image2D w, wp, T;
  int x,y;
  double sigma=sigmin, Tlambda,TGw;

  InitImage2D(&w, in->width, in->height);
  InitImage2D(&wp, in->width, in->height);
  InitImage2D(&T, in->width, in->height);

  setImage2D(&T,2*BIG);         /* set threshold map to out of range value */
  SquareSobel(in,&w);    /* compute gradient */

  eta*=sqrt(5.0)/4.0;
  clipLowImage2D(&w,&w,lambda*lambda*eta*eta); 

  /* remove low edge values, edge now holds weights w(x,y) */

  mulImage2D(&w,in,&wp);    /* Multiply with input to obtain product
                                    image w.p(x,y) */

  Tlambda = sumGreyImage2D(&wp)/sumGreyImage2D(&w);

  recGauss2DFilter( wp.data,     /* Convolve w.p(x,y) with Gaussian */
		    wp.data, 
		    wp.width, 
		    wp.height, 
		    sigma );

  recGauss2DFilter( w.data,    /* Convolve w(x,y) with Gaussian */
		    w.data, 
		    w.width, 
		    w.height, 
		    sigma );
  
  
  for (sigma = sigmin; sigma <sigmax; sigma *= 2.0 )
    { 
      TGw = T_Gw(lambda,sigma,eta);
      //printf("%f\n", TGw);
      for (y=0;y<w.height;y++)
	for (x=0;x<w.width;x++)
	  { 
	    if (T.data[y][x]>BIG)
	      if (w.data[y][x]>TGw)
		T.data[y][x]= wp.data[y][x]/w.data[y][x];
	  }
      recGauss2DFilter( wp.data,     /* Convolve w.p(x,y) with Gaussian */
			wp.data, 
			wp.width, 
			wp.height, 
			sqrt(3)*sigma );

      recGauss2DFilter( w.data,    /* Convolve w(x,y) with Gaussian */
			w.data, 
			w.width, 
			w.height, 
			sqrt(3)*sigma );
      
    }
   
  TGw = T_Gw(lambda,sigma,eta);

  for (y=0;y<w.height;y++)
    for (x=0;x<w.width;x++)
      { 
	if (T.data[y][x]>BIG)
	  if (w.data[y][x]>TGw)
	    T.data[y][x]= wp.data[y][x]/w.data[y][x];
	  else
	    T.data[y][x]= Tlambda;
      }
 
  
  threshImage2D(in,&T,out);

  ExitImage2D(&w);
  ExitImage2D(&wp);
  ExitImage2D(&T);
}

float compareImage2D(Image2D *in1, Image2D *in2) 
{
  int x,y;
  int sum=0;

  for (y=0;y<in1->height;y++)
    {
      for (x=0;x<in1->width;x++)
	sum+= in1->data[y][x] != in2->data[y][x]; 
    }
  return (float)sum/((float)in1->height*in1->width);
}

int main (int argc, char **argv)
{
  Image2D in, noise, ideal, out;
  double sigma=4,lambda=5,eta=3, sigmaMax = 64, sum, sum2, error, slope = 0;
  int idum=5436,x,y, i;
  struct tms timestruct;
  clock_t start;
  long ticksize=sysconf(_SC_CLK_TCK);
  float museconds;
  unsigned long bpp;

  bpp = ReadTIFF(argv[1], &in );
  //  ReadPGM(argv[2], &ideal );

  InitImage2D(&out, in.width, in.height);
  // InitImage2D(&noise, in.width, in.height);


  if (argc>2) lambda=atof(argv[2]); 
  /* read threshold on (square) edge strength */
  if (argc>3) eta=atof(argv[3]);
  if (argc>4) sigma=atof(argv[4]);
  if (argc>5) sigmaMax=atof(argv[5]);

  slope /= (double)in.width;


  times(&timestruct);
  start=timestruct.tms_utime;

  GM_RATS(&in,&out,lambda,eta,sigma,sigmaMax);
	

  times(&timestruct);
  museconds=1E6*(timestruct.tms_utime-start)/(float) (ticksize);
  printf ("%10.4f\n", museconds/1000000);
  WriteTIFF ("gmsrats.tif", &out, 8 );

  ExitImage2D(&in);

  ExitImage2D(&out);

  return;
}
















