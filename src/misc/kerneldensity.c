#include <malloc.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/types.h>
#include <sys/times.h>
#include <unistd.h>
#include <ctype.h>
//#include "avs_io.h"
#include <string.h>
#include <fitsio.h>
#include <float.h> 
#include "kerneldensity.h"
#include "eispack.c"


#define PI 3.141592654


typedef struct { double (*Func)( double sqr_dist);
                 double normalization,
                       range;
               } Kernel;
#define MAX(a,b)  ((a>=b) ? (a) : (b))
#define MIN(a,b)  ((a<=b) ? (a) : (b))


void write_fits_file_basic(char *fnameout,  int scaling,  map3d *map){
  /* +++++++++++++++++++++++++++ */
  /*     FITS Write Function     */
  /* +++++++++++++++++++++++++++ */
  
  fitsfile 	*outfptr, *infptr; 	/* FITS file pointers */
  int 		status = 0; 		/* CFITSIO status value MUST be initialized to zero! */
  int 		type;
  long 		naxes[3]      = {map->x_bins, map->y_bins, map->z_bins};      
  long 		counts[3]     = {map->x_bins, map->y_bins, map->z_bins};      
  long 		offsets[3]    = {1,1,1};
  char 		str1[100];
  strcpy(str1, "!");   		/* '!' symbol makes the output image, if existing, to be overwritten */
  fits_create_file(&outfptr, strcat(str1, fnameout), &status);
  float *buf = malloc(map->x_bins*map->y_bins*map->z_bins*sizeof(float));

  
  type = TFLOAT;
  fits_create_img(outfptr, DOUBLE_IMG, 3, naxes, &status);

  ulong i=0;
  for (ulong z=0; z<map->z_bins; z++)
    for (ulong y=0; y<map->y_bins; y++)
      for (ulong x=0; x<map->x_bins; x++, i++)
	buf[i]=(float) (((double)scaling*(map->map[z][y][x]))/map->data_max);
      
  fits_write_subset(outfptr, type, offsets, counts, buf, &status);
 
  fits_close_file(outfptr, &status);

} /* write_fits_file */


void init_map2d ( map2d *map,
                  double x_min,
		  double x_max,
                  int   x_bins,
                  double y_min,
		  double y_max,
                  int   y_bins)
{ short i;
  map->x_min=x_min;
  map->x_max=x_max;
  map->x_bins=x_bins;
  map->y_min=y_min;
  map->y_max=y_max;
  map->y_bins=y_bins;
  map->data_max=0.0;
  map->data_min=0.0;

  map->map=(double **) calloc(y_bins,sizeof(double *));
  for (i=0;i<y_bins;i++)
    map->map[i]=(double *) calloc(x_bins,sizeof(double));
}


void init_map3d ( map3d *map,
                  double x_min,
		  double x_max,
                  int   x_bins,
                  double y_min,
		  double y_max,
                  int   y_bins,
                  double z_min,
		  double z_max,
                  int   z_bins)
{ short i,j;
  map->x_min=x_min;
  map->x_max=x_max;
  map->x_bins=x_bins;

  map->y_min=y_min;
  map->y_max=y_max;
  map->y_bins=y_bins;

  map->z_min=z_min;
  map->z_max=z_max;
  map->z_bins=z_bins;

  printf("x_min=%f,x_max=%f, \n",map->x_min,map->x_max);
  printf("y_min=%f,y_max=%f, \n",map->y_min,map->y_max);
  printf("z_min=%f,z_max=%f, \n",map->z_min,map->z_max);

  map->data_max=0.0;
  map->data_min=0.0;

  map->map=(double ***) calloc(z_bins,sizeof(double **));
  for (j=0;j<z_bins;j++){
    map->map[j]=(double **) calloc(y_bins,sizeof(double *));
    
    for (i=0;i<y_bins;i++)
      map->map[j][i]=(double *) calloc(x_bins,sizeof(double));
  }
}


void exit_map2d ( map2d *map )
{ short i;
  for (i=0;i<map->y_bins;i++)
    free(map->map[i]);
  free(map->map);
}

void exit_map3d ( map3d *map )
{ short i,j;
  for (i=0;i<map->z_bins;i++){
    for (j=0;j<map->y_bins;j++){
      free(map->map[i][j]);
    }
    free(map->map[i]);
  }
  free(map->map);
}

void exit_regress_map2d ( regress_map2d *map )
{ short i;
  for (i=0;i<map->y_bins;i++)
    { free(map->dens_map[i]);
      free(map->z_map[i]);
    }
  free(map->dens_map);
  free(map->z_map);
}

void exit_map1d ( map1d *map )
{
  free(map->map);
}

void exit_regress_map1d ( regress_map1d *map )
{
  free(map->x_map);
  free(map->y_map);
  free(map->sig_map);
}




void add_one_point_epan ( double       x_val,
                          double       y_val,
                          map2d       *d,
                          double       hx,
                          double       hy    )
/* 
   Add one point using Epanechnikov kernel, without normalization 
   Used for fixed kernel density estimates, which can be post-normalized
   for added speed
*/
 
{ short x,y,x_start;
  double x0,
        cur_y,
        x_step1=(d->x_max- d->x_min)/(hx*(d->x_bins-1)),
        y_step=(d->y_max-
                d->y_min)/(hy*(d->y_bins-1)),
        x_step2=x_step1*x_step1,
        cx_xs,
        dist;
  cur_y=( (d->y_min)-y_val)/hy;
  if (cur_y<(-1))
    { y=(int)((double)(fabs(cur_y)-1)/y_step);
      cur_y+=y*y_step;
    }
  else
    y=0;
  x0=( (d->x_min)-x_val)/hx;
  if (x0<(-1))
    { x_start=(int)((double)(fabs(x0)-1)/x_step1);
      x0+=x_start*x_step1;
    }
  else
    x_start=0;
  x_step1*=x0;
  while (y<d->y_bins && cur_y<=1)
    {
      cx_xs=x_step1;
      x=x_start;
      dist=cur_y*cur_y+x0*x0;
      while (cx_xs<0 && dist>1)
        { dist+=(cx_xs+cx_xs+x_step2);
          cx_xs+=x_step2;
          x++;
        }
      while ( x<d->x_bins && dist<=1)
        { (d->map[y][x])+= (1-dist);   /* no normalization for h */
          dist+=(cx_xs+cx_xs+x_step2);
          cx_xs+=x_step2;
          x++;
        }
      cur_y+=y_step;
      y++;
    }
}


void add_many_points_epan ( double       x_val,
			    double       y_val,
			    double       count,
			    map2d       *d,
			    double       hx,
			    double       hy    )
/* 
   Add multiple points using Epanechnikov kernel at same position, 
   without normalization 
   Used for fixed kernel density estimates, which can be post-normalized
   for added speed
*/
 
{ short x,y,x_start;
  double x0,
        cur_y,
        x_step1=(d->x_max- d->x_min)/(hx*(d->x_bins-1)),
        y_step=(d->y_max-
                d->y_min)/(hy*(d->y_bins-1)),
        x_step2=x_step1*x_step1,
        cx_xs,
        dist;
  cur_y=( (d->y_min)-y_val)/hy;
  if (cur_y<(-1))
    { y=(int)((double)(fabs(cur_y)-1)/y_step);
      cur_y+=y*y_step;
    }
  else
    y=0;
  x0=( (d->x_min)-x_val)/hx;
  if (x0<(-1))
    { x_start=(int)((double)(fabs(x0)-1)/x_step1);
      x0+=x_start*x_step1;
    }
  else
    x_start=0;
  x_step1*=x0;
  while (y<d->y_bins && cur_y<=1)
    {
      cx_xs=x_step1;
      x=x_start;
      dist=cur_y*cur_y+x0*x0;
      while (cx_xs<0 && dist>1)
        { dist+=(cx_xs+cx_xs+x_step2);
          cx_xs+=x_step2;
          x++;
        }
      while ( x<d->x_bins && dist<=1)
        { (d->map[y][x])+= count*(1-dist);   
                                        /* no normalization for h */
          dist+=(cx_xs+cx_xs+x_step2);
          cx_xs+=x_step2;
          x++;
        }
      cur_y+=y_step;
      y++;
    }
}

void add_one_point_epan_3d ( double       x_val,
			     double       y_val,
			     double       z_val,
			     map3d       *d,
			     double       hx,
			     double       hy,
			     double       hz )
/* 
   Add one point using Epanechnikov kernel, without normalization 
   Used for fixed kernel density estimates, which can be post-normalized
*/
{ long x,y,z,x_start,y_start;
  double 
    x0,y0, cur_y, cur_z,
    x_step1=(d->x_max - d->x_min)/(hx*(d->x_bins-1)),
    y_step=(d->y_max - d->y_min)/(hy*(d->y_bins-1)),
    z_step=(d->z_max - d->z_min)/(hz*(d->z_bins-1)),
    x_step2=x_step1*x_step1,
    cx_xs,
    dist;

  // printf("ADD ONE POINT: %lf, %lf, %lf. Width: %lf, %lf, %lf\n", x_val, y_val, z_val, hx, hy, hz);
  // printf("Steps %lf, %lf, %lf.\n", x_step1, y_step, z_step);

  cur_z=( (d->z_min)-z_val)/hz;
  // printf("cur_z %lf, z_min %lf, z_val %lf, z width %lf\n", cur_z, d->z_min, z_val, hz);

  if (cur_z<(-1)){
    z=(long)((double)(fabs(cur_z)-1)/z_step);
    cur_z+=z*z_step;
  }
  else
    z=0;

  y0=( (d->y_min)-y_val)/hy;
  // printf("y0  %lf, y_min %lf, y_val %lf, y width %lf\n", y0, d->y_min, y_val, hy);
  //
  if (y0<(-1)){
    y_start=(long)((double)(fabs(y0)-1)/y_step);
    //   printf("%lf\n", (double)(fabs(y0)-1)/y_step);

    y0+=y_start*y_step;
  }
  else
    y_start=0;

  x0=( (d->x_min)-x_val)/hx;
  //  printf("x0 %lf, x_min %lf, x_val %lf, x width %lf\n", x0, d->x_min, x_val, hx);

  if (x0<(-1)){
    x_start=(long)((double)(fabs(x0)-1)/x_step1);
    //   printf("%lf\n", (double)(fabs(x0)-1)/x_step1);

    x0+=x_start*x_step1;
  }
  else
    x_start=0;

  // printf("Start: x_start: %d, x %lf, x step %lf\n ",x_start, x0, x_step1);
  // printf("Start: y_start: %d, y %lf, y step %lf\n ",y_start, y0, y_step);
  //printf("Start: z_start: %d, z %lf, z step %lf\n ",z, cur_z, z_step);
  x_step1*=x0;

  //
    // printf("Start: new x step %lf\n ", x_step1);

  while (z<d->z_bins && cur_z<=1){
    cur_y = y0;
    y=y_start;
    //  printf("Next z: %d (curr %lf)\n",z, cur_z);

    while (y<d->y_bins && cur_y<=1){
      cx_xs=x_step1;
      x=x_start;
      //  printf("Start x: %d (cx %lf)\n",x, cx_xs);
      //  printf("Next y: %d (curr %lf)\n",y, cur_y);

      dist=cur_z*cur_z+cur_y*cur_y+x0*x0;
      //  printf("Dist: %lf \n",dist);

      while (cx_xs<0 && dist>1){
	dist+=(cx_xs+cx_xs+x_step2);
	//	printf("New Dist: %lf \n",dist);
	cx_xs+=x_step2;
	x++;
	//	printf("New x  %d (cx %lf) \n",x, cx_xs);
      }
      while ( x<d->x_bins && dist<=1){

	(d->map[z][y][x])+= (1-dist);   /* no normalization for h */
	//	printf("New Dmap at x %d, y %d, z %d: %lf \n",x,y,z,1-dist);

	dist+=(cx_xs+cx_xs+x_step2);
	cx_xs+=x_step2;
	x++;
      }
      cur_y+=y_step;
      y++;
    }
    cur_z+=z_step;
    z++;
  }
}

void add_one_point_epan2 ( double       x_val,
                           double       y_val,
                           map2d       *d,
                           double       hx,
                           double       hy    )
{ short x,y,x_start;
  double x0,
        h2=hx*hy,
        cur_y,
        x_step1=(d->x_max-
                d->x_min)/(hx*(d->x_bins-1)),
        y_step=(d->y_max-
                d->y_min)/(hy*(d->y_bins-1)),
        x_step2=x_step1*x_step1,
        cx_xs,
        dist;
  cur_y=( (d->y_min)-y_val)/hy;
  if (cur_y<(-1))
    { y=(int)((double)(fabs(cur_y)-1)/y_step);
      cur_y+=y*y_step;
    }
  else
    y=0;
  x0=( (d->x_min)-x_val)/hx;
  if (x0<(-1))
    { x_start=(int)((double)(fabs(x0)-1)/x_step1);
      x0+=x_start*x_step1;
    }
  else
    x_start=0;
  x_step1*=x0;
  while (y<d->y_bins && cur_y<=1)
    {
      cx_xs=x_step1;
      x=x_start;
      dist=cur_y*cur_y+x0*x0;
      while (cx_xs<0 && dist>1)
        { dist+=(cx_xs+cx_xs+x_step2);
          cx_xs+=x_step2;
          x++;
        }
      while ( x<d->x_bins && dist<=1)
        { (d->map[y][x])+= ((1-dist)/h2); /* normalized for h */
          dist+=(cx_xs+cx_xs+x_step2);
          cx_xs+=x_step2;
          x++;
        }
      cur_y+=y_step;
      y++;
    }
}

void add_many_points_epan2 ( double       x_val,
			     double       y_val,
			     double       count,        
			     map2d       *d,
			     double       hx,
			     double       hy    )
{ short x,y,x_start;
  double x0,
        h2=hx*hy,
        cur_y,
        x_step1=(d->x_max-
                d->x_min)/(hx*(d->x_bins-1)),
        y_step=(d->y_max-
                d->y_min)/(hy*(d->y_bins-1)),
        x_step2=x_step1*x_step1,
        cx_xs,
        dist;
  cur_y=( (d->y_min)-y_val)/hy;
  if (cur_y<(-1))
    { y=(int)((double)(fabs(cur_y)-1)/y_step);
      cur_y+=y*y_step;
    }
  else
    y=0;
  x0=( (d->x_min)-x_val)/hx;
  if (x0<(-1))
    { x_start=(int)((double)(fabs(x0)-1)/x_step1);
      x0+=x_start*x_step1;
    }
  else
    x_start=0;
  x_step1*=x0;
  while (y<d->y_bins && cur_y<=1)
    {
      cx_xs=x_step1;
      x=x_start;
      dist=cur_y*cur_y+x0*x0;
      while (cx_xs<0 && dist>1)
        { dist+=(cx_xs+cx_xs+x_step2);
          cx_xs+=x_step2;
          x++;
        }
      while ( x<d->x_bins && dist<=1)
        { (d->map[y][x])+= count*((1-dist)/h2); /* normalized for h */
          dist+=(cx_xs+cx_xs+x_step2);
          cx_xs+=x_step2;
          x++;
        }
      cur_y+=y_step;
      y++;
    }
}



void add_one_point_epan2_3d ( double       x_val,
			      double       y_val,
			      double       z_val,
			      map3d       *d,
			      double       hx,
			      double       hy,
			      double       hz )
/* 
   Add one point using Epanechnikov kernel, with partial normalization
   (hx*hy*hz)    
   Used for adaptive kernel density estimates, which can be post-normalized
   for number of data and kernel shape, but not for bandwidth.
*/
{ short x,y,z,x_start,y_start;
  double 
    h3=hx*hy*hz,
    x0,y0, cur_y, cur_z,
    x_step1=(d->x_max - d->x_min)/(hx*(d->x_bins-1)),
    y_step=(d->y_max - d->y_min)/(hy*(d->y_bins-1)),
    z_step=(d->z_max - d->z_min)/(hz*(d->z_bins-1)),
    x_step2=x_step1*x_step1,
    cx_xs,
    dist;
  
  cur_z=( (d->z_min)-z_val)/hz;
  if (cur_z<(-1)){
    z=(int)((double)(fabs(cur_z)-1)/z_step);
    cur_z+=z*z_step;
  }
  else
    z=0;

  y0=( (d->y_min)-y_val)/hy;
  if (y0<(-1)){
    y_start=(int)((double)(fabs(y0)-1)/y_step);
    y0+=y_start*y_step;
  }
  else
    y_start=0;

  x0=( (d->x_min)-x_val)/hx;
  if (x0<(-1)){
    x_start=(int)((double)(fabs(x0)-1)/x_step1);
    x0+=x_start*x_step1;
  }
  else
    x_start=0;

  x_step1*=x0;

  while (z<d->z_bins && cur_z<=1){
    cur_y = y0;
    y=y_start;
    while (y<d->y_bins && cur_y<=1){
      cx_xs=x_step1;
      x=x_start;
      dist=cur_z*cur_z+cur_y*cur_y+x0*x0;
      while (cx_xs<0 && dist>1){
	dist+=(cx_xs+cx_xs+x_step2);
	cx_xs+=x_step2;
	x++;
      }
      while ( x<d->x_bins && dist<=1){
	(d->map[z][y][x]) += (1-dist)*(1-dist); // exp(-0.5*dist*dist);   /* normalization for h */
	dist+=(cx_xs+cx_xs+x_step2);
	cx_xs+=x_step2;
	x++;
      }
      cur_y+=y_step;
      y++;
    }
    cur_z+=z_step;
    z++;
  }
}

void init_contours ( contour_rec *contours,
                     int         num_contours,
                     double       min,
                     double       max,
                     double       log_radix,
                     int         min_col,
                     int         max_col       )
{ int i;
  contours->num_contours=num_contours;
  contours->levels=calloc(num_contours,sizeof(contours->levels[0]));
  contours->colours=calloc(num_contours,sizeof(contours->colours[0]));
  for (i=0;i<num_contours;i++)
    { contours->colours[i]=min_col+i*(max_col-min_col)/(num_contours-1);
      contours->levels[i]=min+(double)i*(max-min)/(double)(num_contours-1);
    }
  if ( log_radix!=0 )
    { contours->levels[num_contours-1]=max;
      for (i=num_contours-2; i>=0; i--)
        contours->levels[i]=contours->levels[i+1]/log_radix;
    }
}

void exit_contours ( contour_rec *contours )
{ free(contours->levels);
  free(contours->colours);
}



		      
double prob_from_map2d(  map2d *map,
                         double x,
                         double y,
                         double xstep,
                         double ystep )       
{
  int xbin=(int)((x - map->x_min)/xstep);
  double xratio=(x-map->x_min-(xstep*xbin))/xstep;
  int ybin=(int)((y-map->y_min)/ystep);
  double yratio=(y-map->y_min-(ystep*ybin))/ystep;
  double dens1,dens2;

  dens1 = map->map[ybin][xbin]+( map->map[ybin][xbin+1]
				 -map->map[ybin][xbin]  )*xratio; 
  dens2 = map->map[ybin+1][xbin]+( map->map[ybin+1][xbin+1]
				   -map->map[ybin+1][xbin]  )*xratio;
  return dens1 + ( dens2 - dens1 ) * yratio;
}
double prob_from_map3d(  map3d *map,
                         double x,
                         double y,

                         double z,
                         double xstep,
                         double ystep,
                         double zstep )       
{
 
  int xbin=(int)((x-map->x_min)/xstep);
  double xratio=(x-map->x_min-(xstep*xbin))/xstep;
  int ybin=(int)((y-map->y_min)/ystep);
  double yratio=(y-map->y_min-(ystep*ybin))/ystep;
  int zbin=(int)((z-map->z_min)/zstep);
  double zratio=(z-map->z_min-(zstep*zbin))/zstep;

  double dens1, dens2, dens3, dens4;

  dens1 = map->map[zbin][ybin][xbin]+( map->map[zbin][ybin][xbin+1]
				 -map->map[zbin][ybin][xbin]  )*xratio; 
  dens2 = map->map[zbin][ybin+1][xbin]+( map->map[zbin][ybin+1][xbin+1]
				   -map->map[zbin][ybin+1][xbin]  )*xratio;


  dens3 = dens1 + ( dens2 - dens1 ) * yratio;

  dens1 = map->map[zbin+1][ybin][xbin]+( map->map[zbin+1][ybin][xbin+1]
				 -map->map[zbin+1][ybin][xbin]  )*xratio; 
  dens2 = map->map[zbin+1][ybin+1][xbin]+( map->map[zbin+1][ybin+1][xbin+1]
				   -map->map[zbin+1][ybin+1][xbin]  )*xratio;


  dens4 = dens1 + ( dens2 - dens1 ) * yratio;
  //  printf("%lf \n", x);

  return dens3 + ( dens4 - dens3 ) * zratio;
}



double comp_data_probs_2d ( map2d  *density,
                            double xwidth,
                            double ywidth,
                            int    num_data,
			    double *xdata,
			    double *ydata,
			    double *prob)
{
  int i,k,l;
  double normalization = 1/( (double)num_data*xwidth*ywidth*PI/2.0 ), 
         gmean = 0.0;
  double xstep = (density->x_max-density->x_min)/(double)(density->x_bins - 1);
  double ystep = (density->y_max-density->y_min)/(double)(density->y_bins - 1);


  for (i=0; i<num_data; i++)
    add_one_point_epan(xdata[i],ydata[i],density,xwidth,ywidth);
  density->data_max  =density->map[0][0]*normalization;
  density->data_min = density->map[0][0]*normalization;
  for (k=0; k<density->y_bins; k++)
    for (l=0; l<density->x_bins; l++)
      {
	density->map[k][l] = density->map[k][l]*normalization;
	if (density->map[k][l]>density->data_max) 
	  density->data_max=density->map[k][l];
	else
	  if (density->map[k][l]< density->data_min)
	    density->data_min = density->map[k][l];
      }

  for (i=0; i<num_data; i++)
    {
      prob[i]=prob_from_map2d(density,xdata[i],ydata[i],xstep,ystep);
      gmean += log(prob[i]);
    }

  for (k=0; k<density->y_bins; k++)
    for (l=0; l<density->x_bins; l++)
      density->map[k][l]=0;

  return exp(gmean/num_data);
}


void map3d_normalize( map3d *density, double normalization ){
  int x,y,z;

  density->data_max  =density->map[0][0][0]*normalization;
  density->data_min = density->map[0][0][0]*normalization;
  for (z=0; z<density->z_bins; z++)
    for (y=0; y<density->y_bins; y++)
      for (x=0; x<density->x_bins; x++){
	density->map[z][y][x] = density->map[z][y][x]*normalization;
	if (density->map[z][y][x]>density->data_max) 
	  density->data_max=density->map[z][y][x];
	else
	  if (density->map[z][y][x]< density->data_min)
	    density->data_min = density->map[z][y][x];
      }

}



double *hessian(double ***gvals, ulong *dims, ulong x, ulong y, ulong z){

  double hxx, hyy, hzz, hxy, hxz, hyz;
  // ulong x = p %dims[0];
  //ulong y = (p %(dims[0]*dims[1]))/dims[0];
  //ulong z = (p /dims[0]*dims[1]);
  ulong p = z*dims[0]*dims[1]+y*dims[0]+x;
  ulong xbef = x - 1 >= 0? x-1: dims[0]-1;
  ulong xaf  = x + 1 < dims[0]? x+1: 0;
  ulong ybef = y - 1 >= 0? y-1: dims[1]-1;
  ulong yaf = y + 1 < dims[1]?  y+1: 0;
  ulong zbef = z - 1 >= 0? z-1: dims[2]-1;
  ulong zaf = z + 1 < dims[2]?  z+1: 0;
  static double hessian_arr[9];// = {hxx, hxy, hxz, hxy, hyy, hyz, hxz, hyz, hzz};

  hessian_arr[0] = (double) (gvals[z][y][xbef] + gvals[z][y][xaf] - 2*gvals[z][y][x]);
  // printf("%1.20lf %1.20lf %1.20lf\n",gvals[z][y][xbef], gvals[z][y][xaf], 2*gvals[z][y][x]);

  hessian_arr[4] = (double) (gvals[z][ybef][x] + gvals[z][yaf][x] - 2*gvals[z][y][x]);
  hessian_arr[8] = (double) (gvals[zbef][y][x] + gvals[zaf][y][x] - 2*gvals[z][y][x]);

  hessian_arr[1] = hessian_arr[3] = (double) (gvals[z][ybef][xbef] + gvals[z][yaf][xaf] - gvals[z][ybef][xaf] -gvals[z][yaf][xbef]);
  hessian_arr[2] = hessian_arr[6] = (double) (gvals[zbef][y][xbef] + gvals[zaf][y][xaf] - gvals[zaf][y][xbef] -gvals[zbef][y][xaf]);
  hessian_arr[5] = hessian_arr[7] = (double) (gvals[zbef][ybef][x] + gvals[zaf][yaf][x] - gvals[zbef][yaf][x] -gvals[zaf][ybef][x]);
  //printf("%1.20lf %1.20lf %1.20lf\n", hessian_arr[0], hessian_arr[1], hessian_arr[2]);

  return hessian_arr;
}


double comp_data_probs_3d ( map3d  *density,
                            double xwidth,
                            double ywidth,
                            double zwidth,
                            int    num_data,
			    double *xdata,
			    double *ydata,
			    double *zdata,
			    double *weights,
			    double *prob,
			    double *xshear,
			    double *yshear,
			    double *zshear)
{
  // xwidth *= 0.5;
  //ywidth *= 1.5;
  
  double normalization = 1/sqrt(2*PI); //15.0/( (double)num_data*xwidth*ywidth*zwidth*PI*8.0), 
  double  gmean = 0.0;
  double xstep = (density->x_max-density->x_min)/(double)(density->x_bins - 1);
  double ystep = (density->y_max-density->y_min)/(double)(density->y_bins - 1);
  double zstep = (density->z_max-density->z_min)/(double)(density->z_bins - 1);



  for (unsigned long i=0; i<num_data; i++){
    add_one_point_epan_3d(xdata[i],ydata[i],zdata[i], density,weights[i],weights[i],weights[i]);
  }

  map3d_normalize( density, normalization );
  ulong dims[3] = {density->x_bins,density->x_bins,density->x_bins};
  for (unsigned long i=0; i<num_data; i++)
    {
      prob[i]=prob_from_map3d(density,
			      xdata[i],ydata[i],zdata[i],
			      xstep,ystep,zstep);
      /*   double *hessian_arr = hessian(density->map, dims, xdata[i],ydata[i],zdata[i]);
      //     printf("\n%.20lf %.20lf %.20lf\n", hessian_arr[0], hessian_arr[1], hessian_arr[2]);
      //    printf("%.20lf %.20lf %.20lf\n", hessian_arr[3], hessian_arr[4], hessian_arr[5]);
  //    printf("%.20lf %.20lf %.20lf\n\n", hessian_arr[6], hessian_arr[7], hessian_arr[8]);

      double eigvec[9] = {0};
      double eigval[3] = {0};
      int ierr = rs (3, hessian_arr, eigval, 1, eigvec);
      // printf("%.20lf %.20lf %.20lf\n", eigval[0], eigval[1], eigval[2]);
      xshear[i] = fabs(eigvec[6]);
      yshear[i] = fabs(eigvec[7]);
      zshear[i] = fabs(eigvec[8]);
      */
      gmean += log(prob[i]);
    }

  
  /* printf("gmean = %e \n",exp(gmean/num_data));
     write_fits_file_basic("density_inter.fits", 1, density);*/

  for (unsigned long z=0; z<density->z_bins; z++)
    for (unsigned long y=0; y<density->y_bins; y++)
      for (unsigned long x=0; x<density->x_bins; x++)
      density->map[z][y][x]=0;

  return exp(gmean/num_data);
}



void comp_density_2d ( map2d  *density,
		       double xwidth,
		       double ywidth,
		       double gmean,
		       int    num_data,
		       double *xdata,
		       double *ydata,
		       double *prob)
{
  int i,k,l;
  double normalization = 1/( (double)num_data*PI/2.0 );

  for (i=0; i<num_data; i++)
    add_one_point_epan2(xdata[i],ydata[i],
			density,
			xwidth/sqrt(prob[i]/gmean),
			ywidth/sqrt(prob[i]/gmean));
  density->data_max  =density->map[0][0]*normalization;
  density->data_min = density->map[0][0]*normalization;
  for (k=0; k<density->y_bins; k++)
    for (l=0; l<density->x_bins; l++)
      {
	density->map[k][l] = density->map[k][l]*normalization;
	if (density->map[k][l]>density->data_max) 
	  density->data_max=density->map[k][l];
	else
	  if (density->map[k][l]< density->data_min)
	    density->data_min = density->map[k][l];
      }
}

void comp_density_3d ( map3d  *density,
		       double xwidth,
		       double ywidth,
		       double zwidth,
		       double gmean,
		       int    num_data,
		       double *xdata,
		       double *ydata,
		       double *zdata,
		       double *prob,
		       double *weights,
		       double *dist,
		       ulong *id,
		       double lambda)
{
  int i;
  double normalization =15.0/16.0;//1/sqrt(2*3.14);// * 15.0/( (double)num_data*PI*8.0 );
  bool *visited = calloc(density->z_bins*density->z_bins*density->z_bins, sizeof(bool));
  ulong index = 0;
  ulong dim = density->z_bins;
  printf("%ld \n", num_data);
  for (i=0; i<num_data; i++){
    //if(i != 0 && prob[i] == prob[i-1]) continue;
    //  double norm =MIN(MIN(xshear[i],yshear[i]),zshear[i]);
    //   printf(" %lf, %lf, %lf, %lf \n",norm, xshear[i]/norm, yshear[i]/norm, zshear[i]/norm);
    /*  double dist_curr = dist[id[i]];
	double xdir = fabs(xdata[i]-xdata[id[i]]);
	double ydir = fabs(ydata[i]-ydata[id[i]]);
	double zdir = fabs(zdata[i]-zdata[id[i]]);*/
    if(weights[i] >lambda) continue;
    double xi = round(xdata[i]);
    double yi = round(ydata[i]);
    double zi = round(zdata[i]);
    index = zi*dim*dim+yi*dim+xi;
    //if(weights[i] > 10) weights[i] = 10;
    if(!visited[index]){
      add_one_point_epan2_3d(xi,yi,zi,
			     density,
			     weights[i]+3,
			     weights[i]+3,
			     weights[i]+3
			     );
      visited[index] = 1;
    }
       
    for (ulong j = i+1; j< num_data; j++){
      double xj = round(xdata[j]);
      double yj = round(ydata[j]);
      double zj = round(zdata[j]); 
      double dist_curr = sqrt(pow(xi-xj,2) + pow(yi-yj,2) + pow(zi-zj,2));
      if(dist_curr < 20 && dist_curr > 5){
	double incx = (xi-xj)/dist_curr;
	double incy = (yi-yj)/dist_curr;
	double incz = (zi-zj)/dist_curr;
	double newx = xi;
	double newy = yi;
	double newz = zi;
	double dist_upd = sqrt(pow(newx-xj,2) + pow(newy-yj,2) + pow(newz-zj,2));
       	double frac = dist_upd/dist_curr;
	double oldx = newx;
	double oldy = newy;
	double oldz = newz;
	double wcur = weights[i];
	while(dist_upd > weights[j]/2){
	  // printf("New: %lf, %lf, %lf \n",newx, newy, newz);
	  newx -= incx;
	  newy -= incy;
	  newz -= incz;
	  index = round(newz)*dim*dim+round(newy)*dim+round(newx);
	  dist_upd = sqrt(pow(newx-xj,2) + pow(newy-yj,2) + pow(newz-zj,2));
	  frac = dist_upd/dist_curr;
	  if(!visited[index]){
	    add_one_point_epan2_3d(round(newx),round(newy),round(newz),
				   density,
				   frac*wcur/2+(1-frac)*weights[j]/2+3,
				   frac*wcur/2+(1-frac)*weights[j]/2+3,
				   frac*wcur/2+(1-frac)*weights[j]/2+3
				   );
	    visited[index] = 1;
	  }
	  wcur = frac*wcur;
	    // dist_upd = sqrt(pow(newx-xj,2) + pow(newy-yj,2) + pow(newz-zj,2));
	  oldx = round(newx);
	  oldy = round(newy);
	  oldz = round(newz);
	}
      }
    } 
  }

  map3d_normalize( density, normalization );
}

void interp_density_3d ( map3d  *density,
		       double xwidth,
		       double ywidth,
		       double zwidth,
		       int    num_data,
		       double *xdata,
		       double *ydata,
		       double *zdata,
		       double *weights,
		       double *dist,
			 ulong *id)
{
  int i;
   double normalization = 15.0/( (double)num_data*PI*8.0 );
  for (i=0; i<num_data; i++){
    ulong j = id[i];
    double dist_curr = dist[i];
    double oldw = weights[i];
    double inc = fabs(xdata[i]-xdata[j])/10;
    double newx = xdata[i] + inc;
    double newy = ydata[i] + inc;
    double newz = zdata[i] + inc;
    double oldx, oldy, oldz;
    double new_dist = sqrt(pow(xdata[i]-newx,2) + pow(ydata[i]-newy,2) + pow(zdata[i]-newz,2));
    printf("Distance %lf, %lf\n", dist[i], new_dist);

    while(new_dist < dist_curr){
      printf("Adding at %lf, %lf, %lf\n", newx, newy, newz);
      double frac = new_dist/dist_curr;
      double weib = frac*oldw + (1-frac)*weights[j];
      add_one_point_epan2_3d(newx,newy,newz, density, weib, weib, weib);
      oldw = weib;
      oldx = newx;
      oldy = newy;
      oldz = newz;      
      newx = newx + inc;
      newy = newy + inc;
      newz = newz + inc;
      dist_curr = new_dist;
          double new_dist = sqrt(pow(oldx-newx,2) + pow(oldy-newy,2) + pow(oldz-newz,2));

    }
  }

   map3d_normalize( density, normalization );

}


void findminmax(int numdata, double *data, double *min, double *max){
  int i;

  *max = *min = data[0];
  for (i=1;i<numdata;i++){
    if (data[i]>(*max))
      *max=data[i];
    else
      if (data[i]<(*min))
	*min=data[i];
  }
}


int ImagePGMBinWrite(map2d *map, char *fname)
{
   FILE *outfile;
   int k,l,i;
   unsigned char *buf = malloc(map->x_bins*map->y_bins*sizeof(unsigned char));

   outfile = fopen(fname, "wb");
   if (outfile==NULL) {
      fprintf (stderr, "Error: Can't write the image: %s !", fname);
      return(0);
   }
   fprintf(outfile, "P5\n%d %d\n255\n", map->x_bins, map->y_bins);
 
   i=0;
   for(l=map->y_bins - 1;l>=0;l--)
     for (k=0;k<map->x_bins;k++,i++)
       buf[i] = (unsigned char) ((255.0*sqrt(map->map[l][k]))/sqrt(map->data_max)); 
   fwrite(buf, 1, (size_t)(map->x_bins*map->y_bins), outfile);
 
   fclose(outfile);
   return(1);
} /* ImagePGMBinWrite */


/*void AVSdens3dwrite(int scaling,
		    map3d *map,
		    char *fname){
  FILE *outfile;
  int x,y,z,i;
  unsigned short *buf = malloc(map->x_bins*map->y_bins*map->z_bins*sizeof(unsigned short));
  avs_header avs_head;

  avs_head.ndim = 3;
  avs_head.dim1 = map->x_bins;
  avs_head.dim2 = map->y_bins;
  avs_head.dim3 = map->y_bins;
  avs_head.min_x = map->x_min;
  avs_head.min_y = map->y_min;
  avs_head.min_z = map->y_min;
  avs_head.max_x = map->x_max;
  avs_head.max_y = map->y_max;
  avs_head.max_z = map->z_max;
  avs_head.filetype = 0;
  avs_head.skip = 0;
  avs_head.nspace = 3;
  avs_head.veclen = 1;
  avs_head.dataname[0] = '\0';
  avs_head.datatype = 2;

  outfile = fopen(fname, "wb");
  if (outfile==NULL) {
    fprintf (stderr, "Error: Can't write the image: %s !", fname);
    return;
  }

  avs_write_header(outfile, &avs_head);
  i=0;
  for (z=0; z<map->z_bins; z++)
    for (y=0; y<map->y_bins; y++)
      for (x=0; x<map->x_bins; x++, i++){
	buf[i]=(unsigned short) (((double)scaling*(map->map[z][y][x]))/map->data_max);
      }
   fwrite(buf, sizeof(unsigned short), (size_t)(map->x_bins*map->y_bins*map->z_bins), outfile);
   fclose(outfile); 
}
*/

int main(int argc, char *argv[]){
  char *infilename=argv[1];
  FILE *infile = fopen(infilename,"r");
  char *skelname;
  FILE *skelfile = NULL;
  ulong lambda;
  if(argc >= 3){
    skelname = argv[3];
    skelfile = fopen(skelname, "r");
    lambda = atoi(argv[4]);
  } else
    lambda = 0;
  char *outfilename = argv[2];
  ulong datapos=0, dataskel=0, numdata,numcols=3, i,j;
  double **data;
  double *minima,*maxima, *prob, *xshear, *yshear, *zshear;
  double gmean,hscale;
  map2d density;
  map3d density3;
  clock_t start;
  struct tms tstruct;
  long tickspersec = sysconf(_SC_CLK_TCK);  
  double musec;
  int n_dims = 3;
  int k=0,l=1,m=2;
  long dims[3];
  
  if(argc>5){
    k=atoi(argv[5])-1; 
  }
  if(argc>6){
    l=atoi(argv[6])-1; 
  }
  if(argc>7){
    m=atoi(argv[7])-1; 
  }

  ulong tempp;
  data = (double **)calloc(numcols, sizeof(double*));
  minima= (double *)calloc(numcols,sizeof(double));
  maxima= (double *)calloc(numcols,sizeof(double));

  printf("Reading position file from disccoman\n");
  fscanf(infile,"%ld\n",&datapos);
  printf("Number of positions in %s: %ld\n", argv[1], datapos);
  if(skelfile != NULL){
    printf("Also reading in the skeleton file pfskel for regions larger than %ld\n", lambda);
    fscanf(skelfile,"%ld\n",&dataskel);
    printf("Number of segments in %s: %ld\n", argv[2], dataskel/2);
    fscanf(skelfile, "0 %ld 0 %ld 0 %ld\n", &tempp, &tempp,&tempp);
  }
  fscanf(infile,"0 %ld 0 %ld 0 %ld\n",dims, dims+1, dims+2);
  printf("Dimension: %ld, %ld, %ld\n", dims[0], dims[1], dims[2]);
  numdata = datapos+dataskel;
  printf("All data to kernelize: %ld\n", numdata);
      
  for (j=0;j<numcols;j++){
    data[j]= (double *)calloc(numdata,sizeof(double));
  }
  // prob = (double *)calloc(numdata,sizeof(double));
  // xshear = (double *)calloc(numdata,sizeof(double));
  // yshear = (double *)calloc(numdata,sizeof(double));
  // zshear = (double *)calloc(numdata,sizeof(double));

  double *weights = (double *)calloc(numdata, sizeof(double));
  double *seg = (double *) calloc(numdata, sizeof(double));
  double temp;
  double *dist_arr = (double *) calloc(numdata, sizeof(double));
  ulong *id = (ulong *) calloc(numdata, sizeof(double));
  for (i = 0; i<numdata; i++){
    dist_arr[i] = -1;
  }
  for (i = 0; i<datapos; i++){
    //for (j=0;j<numcols;j++)
    // fscanf(infile,"%lf %lf %lf %lf\n",data[0]+i,data[1]+i,data[2]+i, weights+i);
    fscanf(infile,"%lf %lf %lf %lf\n",data[0]+i,data[1]+i,data[2]+i, weights+i);
    // printf("%lf %lf %lf %lf\n",*(data[0]+i),*(data[1]+i),*(data[2]+i), *(weights+i));
  }
  if(skelfile != NULL){
    for (i = datapos; i<numdata; i++){
      //for (j=0;j<numcols;j++)
      // fscanf(infile,"%lf %lf %lf %lf\n",data[0]+i,data[1]+i,data[2]+i, weights+i);
      fscanf(skelfile,"%lf %lf %lf %ld\n",data[0]+i,data[1]+i,data[2]+i, seg+i);
      // printf("%lf, %lf %lf\n", data[0][i], data[1][i], data[2][i]);
      weights[i] = (double) pow(lambda,1.0/3.0)/2;
      // printf("%lf %lf %lf %lf\n",*(data[0]+i),*(data[1]+i),*(data[2]+i), *(weights+i));
    }
  }
  printf("Done reading data\n");

  /*  for (i = 0; i<numdata; i++){
    for (j = i+1; j< numdata; j++){
      double dist_curr = sqrt(pow(data[0][i]-data[0][j],2) + pow(data[1][i]-data[1][j],2) + pow(data[2][i]-data[2][j],2));
      //  printf("%d %f \n", i, dist_curr);
      if(dist_arr[i] == -1 || dist_curr < dist_arr[i]){
	dist_arr[i] = dist_curr;
	id[i] = j;
      }
      if(dist_arr[j] == -1 || dist_curr < dist_arr[j]){
	dist_arr[j] = dist_curr;
	id[j] = i;
      }
    }
  }
  */

 
  maxima[k] = maxima[l]= maxima[m] = dims[0]-1;
  if (n_dims==2){
    init_map2d(&density,minima[k]-0.06125*(maxima[k]-minima[k]),
	       maxima[k]+0.06125*(maxima[k]-minima[k]),512,
	       minima[l]-0.06125*(maxima[l]-minima[l]),
	       maxima[l]+0.06125*(maxima[l]-minima[l]),512);
    
    
    start = times(&tstruct);
  
    hscale = 6.0/sqrt(numdata); 
    
    gmean = comp_data_probs_2d(&density,
			       hscale*(maxima[k]-minima[k]),
			       hscale*(maxima[l]-minima[l]),
			       numdata, data[k], data[l], prob);
    
    
    comp_density_2d ( &density,hscale*(maxima[k]-minima[k]),
		      hscale*(maxima[l]-minima[l]),gmean,numdata,
		      data[k],data[l],prob);
    
    musec = (double)(times(&tstruct) - start)/((double)tickspersec);
    
    printf("wall-clock time: %f s\n",musec);
    
    ImagePGMBinWrite(&density,"density.pgm");
  } else {
    /* init_map3d(&density3,minima[k]-0.06125*(maxima[k]-minima[k]),
	       maxima[k]+0.06125*(maxima[k]-minima[k]),dims[0]*1,
	       minima[l]-0.06125*(maxima[l]-minima[l]),
	       maxima[l]+0.06125*(maxima[l]-minima[l]),dims[1]*1,
	       minima[m]-0.06125*(maxima[m]-minima[m]),
	       maxima[m]+0.06125*(maxima[m]-minima[m]),dims[2]*1);*/
   
     init_map3d(&density3,minima[k]-0.0*(maxima[k]-minima[k]),
	       maxima[k]+0.0*(maxima[k]-minima[k]),dims[0],
	       minima[l]-0.0*(maxima[l]-minima[l]),
	       maxima[l]+0.0*(maxima[l]-minima[l]),dims[1],
	       minima[m]-0.0*(maxima[m]-minima[m]),
		maxima[m]+0.0*(maxima[m]-minima[m]),dims[2]);
   
    start = times(&tstruct);
  
    hscale = 2.0/sqrt(numdata); 
  
    /*    gmean = comp_data_probs_3d(&density3,
			       hscale*(maxima[k]-minima[k]),
			       hscale*(maxima[l]-minima[l]),
			       hscale*(maxima[m]-minima[m]),
			       numdata, data[k], data[l], data[m], weights, prob, xshear, yshear, zshear);
    */
    /* for (unsigned long i=0; i<numdata; i++)
      {
	//    prob[i]=10;
	//	gmean += log(prob[i]);
	  printf("prob %lf\t", prob[i]);

	  }/*/
    
    comp_density_3d ( &density3,
		      hscale*(maxima[k]-minima[k]),
		      hscale*(maxima[l]-minima[l]),
		      hscale*(maxima[m]-minima[m]),
		      gmean,numdata,
		      data[k],data[l],data[m],seg, weights, dist_arr, id, lambda != 0 ? (double) pow(lambda,1.0/3.0):pow(2,63));

    /*   interp_density_3d(&density3,  hscale*(maxima[k]-minima[k]),
		      		      hscale*(maxima[l]-minima[l]),
		      hscale*(maxima[m]-minima[m]),
		      numdata,
		      data[k],data[l],data[m], weights, dist_arr, id);*/
    
    musec = (double)(times(&tstruct) - start)/((double)tickspersec);
    
    printf("wall-clock time: %f s\n",musec);
    write_fits_file_basic(outfilename,1, &density3);
    //  AVSdens3dwrite(4095,&density3,"density.fld");
  }
  return 0;
}



