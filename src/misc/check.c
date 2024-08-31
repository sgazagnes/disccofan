#define _GNU_SOURCE
#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <float.h>
#include <string.h>
#include <stdarg.h>
#include <ctype.h>
#include <sys/wait.h>

#include "io.h"

#define RED  "\x1B[91m"
#define GRN  "\x1B[92m"
#define YEL  "\x1B[93m"
#define BLU  "\x1B[94m"
#define MAG  "\x1B[95m"
#define CYN  "\x1B[96m"
#define WHT  "\x1B[97m"
#define RESET "\033[0m"

typedef uint8_t ubyte;
typedef uint16_t ushort;
typedef uint32_t uint;
typedef uint64_t ulong;
typedef int64_t idx;
typedef uint value_t;

//#define FLOAT_TYPE 1
//typedef float value;
enum {
  IO = 0, FILTER, CSL, PATTERN, ALL = INT_MAX
};

bool equals(const char *a, const char *b) {
  return (strcmp(a, b) == 0);
}


FILE *log_open(const char *filename)
{
  FILE *LOG;
  char cmnd[1000];

  sprintf(cmnd, "echo git commit = `git rev-parse HEAD` > %s", filename);
  system(cmnd);

  LOG = fopen(filename, "a");
  if (!LOG){
    return LOG;
  }

  /* Do more stuff in the future */
  
  return LOG;
}

const char *get_filename_ext(const char *filename) {
    const char *dot = strrchr(filename, '.');
    if(!dot || dot == filename) return "";
    return dot + 1;
}

int compare_txt_files(FILE *fp1, FILE *fp2)
{
  char ch1 = getc(fp1);
  char ch2 = getc(fp2);
  int error = 0, pos = 0, line = 1;
  while (ch1 != EOF && ch2 != EOF)
    {
      pos++;

      if (ch1 == '\n' && ch2 == '\n')
        {
	  line++;
	  pos = 0;
        }

      if (ch1 != ch2)
        {
	  error++;
	  printf("Line Number : %d \t Error "
		 "Position : %d \n", line, pos);
        }
      ch1 = getc(fp1);
      ch2 = getc(fp2);
    }

  if(error)
    printf("\x1B[91m Errors spotted \033[0m\n" );
  else
    printf("\x1B[92m No differences spotted \033[0m\n");

  return error;
}


void compare_files(const char* filein, const char* extin, const char* fileout, const char* extout, int n_images, const char *dataset_in, const char *dataset_out){
  double *im_out, *im_in, *im_inb;
  //char *extin = get_filename_ext(filein);
  // char *extout = get_filename_ext(fileout);

  if(!strcmp(extin, "txt")){
    char *filename_in, *filename_out;
    asprintf(&filename_in, "%s.%s", filein, extin);
    asprintf(&filename_out, "%s.%s", fileout, extout);
    FILE *fp1 = fopen(filename_in, "r");
    FILE *fp2 = fopen(filename_out, "r");
    if (fp1 == NULL || fp2 == NULL)
      {
	printf("Error : Files not open\n");
	exit(0);
      }
 
    int error = compare_txt_files(fp1, fp2);
    fclose(fp1);
    fclose(fp2);
  } else{
    ulong dims_A[3] = {1,1,1};
    ulong dims_B[3] = {1,1,1};
    int   bitpix_A, bitpix_B;
    float diff = 0, maxdiff = 0;
    char *prefix;
    for (int i = 0; i < n_images; i++) {
      n_images > 1 ? asprintf(&prefix, "%s-%d", filein, i) :
	asprintf(&prefix, "%s", filein); 
      im_in = read_file(prefix, extin, dataset_in, dims_A, &bitpix_A);	  		
      int n_dimsA = dims_A[2] > 1 ? 3 : 2;
      n_images > 1 ? asprintf(&prefix, "%s-%d", fileout,i) :
	asprintf(&prefix, "%s", fileout); 
      im_inb = read_file(prefix, extout, dataset_out, dims_B, &bitpix_B);	  	      
      int n_dimsB = dims_B[2] > 1 ? 3 : 2;	
      if(n_dimsA != n_dimsB){
	printf("Error: N dimensions not matching in A (%d) and B (%d)",n_dimsA, n_dimsB);
	exit(0);
      }
      if(bitpix_A != bitpix_B)
	printf("Warn: different bit depth in the A (%d) and B (%d)",bitpix_A, bitpix_B);
      if(dims_A[0] != dims_B[0] || dims_A[1] != dims_B[1] || dims_A[2] != dims_B[2] ){
	printf("Error: dimensions not matching in A (W:%lu, H:%lu, D:%lu) and B (W:%lu, H:%lu, D:%lu)", dims_A[0], dims_A[1], dims_A[2], dims_B[0], dims_B[1], dims_B[2]);
	exit(0);
      }

      ulong size = dims_A[0]*dims_A[1]*dims_A[2];
      im_out = malloc(size * sizeof(double));
	
      for(size_t k = 0; k < size; k++){
	im_out[k] = bitpix_A < 0 ? (double) fabs( im_in[k] -  im_inb[k]) : (double) abs( im_in[k] -  im_inb[k]);
	if(im_out[k] > 0){
	  diff = 1;
	  maxdiff = im_out[k] > maxdiff ? im_out[k] : maxdiff;
	}
      }

      // n_images > 1 ? asprintf(&prefix, "Compare-%d", i) :
	//	asprintf(&prefix, "Compare"); 
      // write_file(prefix, "fits", "data", dims_A, im_out, -32, n_dimsA);
    }
    if(diff)
      printf("\x1B[91m Error: differences spotted (Max Diff is %f) \033[0m\n", maxdiff);
    else
      printf("\x1B[92m No differences spotted \033[0m\n");
  }
}


int main(int argc, char** argv) {

  char *c, *mode;
  int choice,i = 0;
  char cmnd[1000];
  FILE *LOG;
  
  argc <= 1 ? asprintf(&mode, "all"): asprintf(&mode,  argv[1]);

  for (char *p = mode; *p != '\0'; ++p) {
    *p = toupper(*p);
  }

  if (equals(mode, "IO")) {
    choice = IO;
  } else if (equals(mode, "FILTER")) {
    choice = FILTER;
  } else if (equals(mode, "CSL")) {
    choice = CSL;
  } else if (equals(mode, "PATTERN")) {
    choice = PATTERN;
  } else if (equals(mode, "ALL")) {
    choice = ALL;
  } else {
    fprintf(stdout, RED "No valid mode level supplied (choices: IO, FILTER, DP, PATTERN, ALL (default))"RESET);
    fflush(stdout);
    exit(0);
  }

  //  LOG = log_open("check_files/check.log");
  
  /* if (!LOG){
     fprintf(stderr, "check.c: Unable to open log file\n Aborting...\n");
     return -1;
     }*/
 
  fprintf(stdout, "CHECK DISCCOFAN VALIDITY, MODE:" GRN "  %s\n"RESET, mode);  fflush(stdout);
  // fprintf(LOG, "CHECK DISCCOMAN VALIDITY, MODE: %s\n", mode);


  if(choice == IO || choice == ALL){
    fprintf(stdout, CYN "\n\n Checking I/O: start\n" RESET);  fflush(stdout);
    fprintf(stdout,  "\n PGM to PGM, single file \n" );  fflush(stdout);
    system("mpirun -np 1 ./disccofan -g 1,1,1 --inprefix check_files/lena512 --intype pgm -l 0 --threads 1 -v off");
    if( access("out.pgm", F_OK ) != -1){
      compare_files("check_files/lena512","pgm", "out","pgm", 1, NULL, NULL);
      //system("./polyImage compare check_files/lena512 pgm out pgm 1");
      system("rm *pgm *fits");
    } else {
      fprintf(stdout, RED "Error, output files not written \n" RESET);  fflush(stdout);
      exit(0);
    }
    
    fprintf(stdout,  "\n PGM to PGM, single to chunks\n" );  fflush(stdout);
    system("mpirun -np 4 ./disccofan -g 2,2,1 --inprefix check_files/lena512 --intype pgm -l 0 --threads 1 --outfile 1 -v off");
    if( access("out-0.pgm", F_OK ) != -1){
      compare_files("check_files/tile","pgm", "out","pgm", 4, NULL, NULL);

      //  system("./polyImage compare check_files/tile pgm out pgm 4");
      system("rm *pgm *fits");
    } else {
      fprintf(stdout, RED "Error, output files not written \n" RESET);  fflush(stdout);
      exit(0);
    }
    
    fprintf(stdout,  "\n PGM to PGM, chunks to chunks\n" );  fflush(stdout);
    system("mpirun -np 4 ./disccofan -g 2,2,1 --inprefix check_files/tileov --intype pgm --infile 1 -l 0 --threads 1 --outfile 1 -v off --overlap 1");
    if( access("out-0.pgm", F_OK ) != -1){
      // system("./polyImage compare check_files/tile pgm out pgm 4");
      compare_files("check_files/tile","pgm", "out","pgm", 4, NULL, NULL);
	 
      system("rm *pgm *fits");
    } else {
      fprintf(stdout, RED "Error, output files not written \n" RESET);  fflush(stdout);
      exit(0);
    }


    
    fprintf(stdout,  "\n\n\n PGM to FITS, single file\n" );  fflush(stdout);
    system("mpirun -np 1 ./disccofan -g 1,1,1 --inprefix check_files/lena512 --intype pgm -l 0 --threads 1 --outtype fits -v off");
    if( access("out.fits", F_OK ) != -1){
      //    system("./polyImage compare check_files/lena512 pgm out fits 1");
      compare_files("check_files/lena512","pgm", "out","fits",1, NULL, NULL);

      system("rm *fits");
    } else {
      fprintf(stdout, RED "Error, output files not written \n" RESET);  fflush(stdout);
      exit(0);
    }
    
    fprintf(stdout,  "\n PGM to FITS, single to chunks\n" );  fflush(stdout);
    system("mpirun -np 4 ./disccofan -g 2,2,1 --inprefix check_files/lena512 --intype pgm -l 0 --threads 1 --outtype fits --outfile 1 -v off");
    if( access("out-0.fits", F_OK ) != -1){
      // system("./polyImage compare check_files/tile pgm out fits 4");
      compare_files("check_files/tile","pgm", "out","fits", 4, NULL, NULL);
	 
      system("rm *fits");
    } else {
      fprintf(stdout, RED "Error, output files not written \n" RESET);  fflush(stdout);
      exit(0);
    }
    
    fprintf(stdout,  "\n PGM to FITS, chunks to single \n" );  fflush(stdout);
    system("mpirun -np 4 ./disccofan -g 2,2,1 --inprefix check_files/tileov --intype pgm --infile 1 -l 0 --threads 1 --outtype fits --outfile 0 -v off --overlap 1");
    if( access("out.fits", F_OK ) != -1){
      // system("./polyImage compare check_files/lena512 pgm out fits 1");
      compare_files("check_files/lena512","pgm", "out","fits", 1, NULL, NULL);

      system("rm *fits");
    } else {
      fprintf(stdout, RED "Error, output files not written \n" RESET);  fflush(stdout);
      exit(0);
    }
    
    fprintf(stdout,  "\n PGM to FITS, chunks to chunks \n" );  fflush(stdout);
    system("mpirun -np 4 ./disccofan -g 2,2,1 --inprefix check_files/tileov --intype pgm --infile 1 -l 0 --threads 1 --outtype fits --outfile 1 -v off --overlap 1");
    if( access("out-0.fits", F_OK ) != -1){
      // system("./polyImage compare check_files/tile pgm out fits 4");
      compare_files("check_files/tile","pgm", "out","fits", 4, NULL, NULL);

      system("rm *fits");
    } else {
      fprintf(stdout, RED "Error, output files not written \n" RESET);  fflush(stdout);
      exit(0);
    } 

    
    fprintf(stdout,  "\n\n\n PGM to H5, single file\n" );  fflush(stdout);
    system("mpirun -np 1 ./disccofan -g 1,1,1 --inprefix check_files/lena512 --intype pgm -l 0 --threads 1 --outtype h5 -v off");
    if( access("out.h5", F_OK ) != -1){
      // system("./polyImage compare check_files/lena512 pgm out h5 1 p filter_Area_0");
      compare_files("check_files/lena512","pgm", "out","h5", 1, NULL, "filter_Area_0");

      system("rm *fits *h5");
    } else {
      fprintf(stdout, RED "Error, output files not written \n" RESET);  fflush(stdout);
      exit(0);
    }
    
    fprintf(stdout,  "\n PGM to H5, single to chunks\n" );  fflush(stdout);
    system("mpirun -np 4 ./disccofan -g 2,2,1 --inprefix check_files/lena512 --intype pgm -l 0 --threads 1 --outtype h5 --outfile 1 -v off");
    if( access("out-0.h5", F_OK ) != -1){
      //      system("./polyImage compare check_files/tile pgm out h5 4 p filter_Area_0");
      compare_files("check_files/tile","pgm", "out","h5", 4, NULL, "filter_Area_0");

      system("rm *fits *h5");
    } else {
      fprintf(stdout, RED "Error, output files not written \n" RESET);  fflush(stdout);
      exit(0);
    }
    
    fprintf(stdout,  "\n PGM to H5, chunks to single\n" );  fflush(stdout);
    system("mpirun -np 4 ./disccofan -g 2,2,1 --inprefix check_files/tileov --intype pgm --infile 1 -l 0 --threads 1 --outtype h5 --outfile 0 -v off");
    if( access("out.h5", F_OK ) != -1){
      //  system("./polyImage compare check_files/lena512 pgm out h5 1 p filter_Area_0");
        compare_files("check_files/lena512","pgm", "out","h5", 1, NULL, "filter_Area_0");
  
      system("rm *fits *h5");
    } else {
      fprintf(stdout, RED "Error, output files not written \n" RESET);  fflush(stdout);
      exit(0);
    }
    
    fprintf(stdout,  "\n PGM to H5, chunks to chunks \n" );  fflush(stdout);
    system("mpirun -np 4 ./disccofan -g 2,2,1 --inprefix check_files/tileov --intype pgm --infile 1 -l 0 --threads 1 --outtype h5 --outfile 1 -v off --overlap 1");
    if( access("out-0.h5", F_OK ) != -1){
      // system("./polyImage compare check_files/tile pgm out h5 4 p filter_Area_0");
        compare_files("check_files/tile","pgm", "out","h5", 4, NULL, "filter_Area_0");

      system("rm *fits *h5");
    } else {
      fprintf(stdout, RED "Error, output files not written \n" RESET);  fflush(stdout);
      exit(0);
    }
      
    fprintf(stdout,  "\n FITS to FITS, single file (2D) \n" );  fflush(stdout);
    system("mpirun -np 1 ./disccofan -g 1,1,1 --inprefix check_files/lena512 --intype fits -l 0 --threads 1 --outtype fits  -v off");
    if( access("out.fits", F_OK ) != -1){
      //system("./polyImage compare check_files/lena512 fits out fits 1 ");
      compare_files("check_files/lena512","pgm", "out","fits", 1, NULL, NULL);

      system("rm *fits");
    } else {
      fprintf(stdout, RED "Error, output files not written \n" RESET);  fflush(stdout);
      exit(0);
    }
    
    fprintf(stdout,  "\n FITS to FITS, single file to chunks (2D) \n" );  fflush(stdout);
    system("mpirun -np 4 ./disccofan -g 2,2,1 --inprefix check_files/lena512 --intype fits -l 0 --threads 1 --outtype fits --outfile 1 -v off");
    if( access("out-0.fits", F_OK ) != -1){
      // system("./polyImage compare check_files/tile pgm out fits 4 ");
      compare_files("check_files/tile","pgm", "out","fits", 4, NULL, NULL);

      system("rm *fits");
    } else {
      fprintf(stdout, RED "Error, output files not written \n" RESET);  fflush(stdout);
      exit(0);
    }
    
    fprintf(stdout,  "\n FITS to FITS, chunks to chunks (2D) \n" );  fflush(stdout);
    system("mpirun -np 4 ./disccofan -g 2,2,1 --inprefix check_files/tileov --intype fits --infile 1 -l 0 --threads 1 --outtype fits --outfile 1 -v off --overlap 1");
    if( access("out-0.fits", F_OK ) != -1){
      //system("./polyImage compare check_files/tile pgm out fits 4 ");
      compare_files("check_files/tile","pgm", "out","fits", 4, NULL, NULL);

      system("rm *fits");
    } else {
      fprintf(stdout, RED "Error, output files not written \n" RESET);  fflush(stdout);
      exit(0);
    }
    
    fprintf(stdout,  "\n FITS to FITS, chunks to single (2D) \n" );  fflush(stdout);
    system("mpirun -np 4 ./disccofan -g 2,2,1 --inprefix check_files/tileov --intype fits --infile 1 -l 0 --threads 1 --outtype fits --outfile 0 -v off --overlap 1");
    if( access("out.fits", F_OK ) != -1){
      // system("./polyImage compare check_files/lena512 fits out fits 1 ");
      compare_files("check_files/lena512","pgm", "out","fits", 1, NULL, NULL);

      system("rm *fits");
    } else {
      fprintf(stdout, RED "Error, output files not written \n" RESET);  fflush(stdout);
      exit(0);
    }


      
    fprintf(stdout,  "\n\n\n FITS to FITS, single file (3D) \n" );  fflush(stdout);
    system("mpirun -np 1 ./disccofan -g 1,1,1 --inprefix check_files/Check16b --intype fits -l 0 --threads 1 --outtype fits  -v off");
    if( access("out.fits", F_OK ) != -1){
      // system("./polyImage compare check_files/Check16b fits out fits 1 ");
      compare_files("check_files/Check16b","fits", "out","fits", 1, NULL, NULL);

      system("rm *fits");
    } else {
      fprintf(stdout, RED "Error, output files not written \n" RESET);  fflush(stdout);
      exit(0);
    }
    
    fprintf(stdout,  "\n FITS to FITS, single file to chunks (3D) \n" );  fflush(stdout);
    system("mpirun  -np 8 ./disccofan -g 2,2,2 --inprefix check_files/Check16b --intype fits -l 0 --threads 1 --outtype fits --outfile 1 -v off");
    if( access("out-0.fits", F_OK ) != -1){
      //system("./polyImage compare check_files/Checktile fits out fits 8 ");
      compare_files("check_files/Checktile","fits", "out","fits", 8, NULL, NULL);

      system("rm *fits");
    } else {
      fprintf(stdout, RED "Error, output files not written \n" RESET);  fflush(stdout);
      exit(0);
    }

    fprintf(stdout,  "\n FITS to FITS, chunks to chunks (3D) \n" );  fflush(stdout);
    system("mpirun  -np 8 ./disccofan -g 2,2,2 --inprefix check_files/Checktileov --intype fits --infile 1 -l 0 --threads 1 --outtype fits --outfile 1 -v off");
    if( access("out-0.fits", F_OK ) != -1){
      //     system("./polyImage compare check_files/Checktile fits out fits 8 ");
      compare_files("check_files/Checktile","fits", "out","fits", 8, NULL, NULL);

      system("rm *fits");
    } else {
      fprintf(stdout, RED "Error, output files not written \n" RESET);  fflush(stdout);
      exit(0);
    }

    fprintf(stdout,  "\n FITS to FITS, chunks to single (3D) \n" );  fflush(stdout);
    system("mpirun -np 8 ./disccofan -g 2,2,2 --inprefix check_files/Checktileov --intype fits --infile 1 -l 0 --threads 1 --outtype fits --outfile 0 -v off --overlap 1");
    if( access("out.fits", F_OK ) != -1){
      //  system("./polyImage compare check_files/Check16b fits out fits 1 ");
      compare_files("check_files/Check16b","fits", "out","fits", 1, NULL, NULL);

      system("rm *fits");
    } else {
      fprintf(stdout, RED "Error, output files not written \n" RESET);  fflush(stdout);
      exit(0);
    }

      
    fprintf(stdout,  "\n\n\n FITS to H5, single file (3D) \n" );  fflush(stdout);
    system("mpirun -np 1 ./disccofan -g 1,1,1 --inprefix check_files/Check16b --intype fits -l 0 --threads 1 --outtype h5  -v off");
    if( access("out.h5", F_OK ) != -1){
      //   system("./polyImage compare check_files/Check16b fits out h5 1 p filter_Area_0 ");
      compare_files("check_files/Check16b","fits", "out", "h5", 1, NULL, "filter_Area_0");

      system("rm *fits *h5");
    } else {
      fprintf(stdout, RED "Error, output files not written \n" RESET);  fflush(stdout);
      exit(0);
    }

    fprintf(stdout,  "\n FITS to H5, single file to chunks (3D) \n" );  fflush(stdout);
    system("mpirun  -np 8 ./disccofan -g 2,2,2 --inprefix check_files/Check16b --intype fits -l 0 --threads 1 --outtype h5 --outfile 1 -v off");
    if( access("out-0.h5", F_OK ) != -1){
      //  system("./polyImage compare check_files/Checktile fits out h5 8 p filter_Area_0");
      compare_files("check_files/Checktile","fits", "out", "h5", 8, NULL, "filter_Area_0");

      system("rm *fits *h5");
    } else {
      fprintf(stdout, RED "Error, output files not written \n" RESET);  fflush(stdout);
      exit(0);
    }

    fprintf(stdout,  "\n FITS to H5, chunks to chunks (3D) \n" );  fflush(stdout);
    system("mpirun  -np 8 ./disccofan -g 2,2,2 --inprefix check_files/Checktileov --intype fits --infile 1 -l 0 --threads 1 --outtype h5 --outfile 1 -v off --overlap 1 ");
    if( access("out-0.h5", F_OK ) != -1){
      // system("./polyImage compare check_files/Checktile fits out h5 8  p filter_Area_0");
      compare_files("check_files/Checktile","fits", "out", "h5", 8, NULL, "filter_Area_0");

      system("rm *fits *h5");
    } else {
      fprintf(stdout, RED "Error, output files not written \n" RESET);  fflush(stdout);
      exit(0);
    }

    fprintf(stdout,  "\n FITS to H5, chunks to single (3D) \n" );  fflush(stdout);
    system("mpirun  -np 8 ./disccofan -g 2,2,2 --inprefix check_files/Checktileov --intype fits --infile 1 -l 0 --threads 1 --outtype h5 --outfile 0 -v off --overlap 1");
    if( access("out.h5", F_OK ) != -1){
      // system("./polyImage compare check_files/Check16b fits out h5 1 p filter_Area_0 ");
      compare_files("check_files/Check16b","fits", "out", "h5", 1, NULL, "filter_Area_0");

      system("rm *fits *h5");
    } else {
      fprintf(stdout, RED "Error, output files not written \n" RESET);  fflush(stdout);
      exit(0);
    }
  }

  if(choice == FILTER || choice == ALL){
    fprintf(stdout, CYN "\n\n Checking Filtering: start\n" RESET);  fflush(stdout);
    fprintf(stdout,  "\n\n\n Lena (512x512 2D 8bit), 1 MPI proc, 1 thread \n" );  fflush(stdout);
    system("mpirun -np 1 ./disccofan -g 1,1,1 --inprefix check_files/lena512 --intype pgm -o test --threads 1 -v off");
    if( access("out.fits", F_OK ) != -1){
      // system("./polyImage compare check_files/lena_check fits out fits 1");
      compare_files("check_files/lena_check","fits", "out", "fits", 1, NULL, NULL);

      system("rm *fits");
    } else {
      fprintf(stdout, RED "Error, output files not written \n" RESET);  fflush(stdout);
      exit(0);
    }
   

    fprintf(stdout,  "\n\n\n Lena (512x512 2D 8bit ), 32 MPI proc, 1 thread \n" );  fflush(stdout);
    system("mpirun  -np 32 ./disccofan -g 8,4,1 --inprefix check_files/lena512 --intype pgm -o test --threads 1 -v off --outtype fits");
    if( access("out.fits", F_OK ) != -1){
      //   system("./polyImage compare check_files/lena_check fits out fits 1");
      compare_files("check_files/lena_check","fits", "out", "fits", 1, NULL, NULL);

      system("rm *fits");
    } else {
      fprintf(stdout, RED "Error, output files not written \n" RESET);  fflush(stdout);
      exit(0);
    }
    
    
    fprintf(stdout,  "\n\n\n Lena (512x512 2D 8bit), 4 MPI proc, 8 thread \n" );  fflush(stdout);
    system("mpirun  -np 4 ./disccofan -g 2,2,1 --inprefix check_files/lena512 --intype pgm -o test --threads 8 -v off --outtype fits");
    if( access("out.fits", F_OK ) != -1){
      //  system("./polyImage compare check_files/lena_check fits out fits 1");
      compare_files("check_files/lena_check","fits", "out", "fits", 1, NULL, NULL);

      system("rm *fits");
    } else {
      fprintf(stdout, RED "Error, output files not written \n" RESET);  fflush(stdout);
      exit(0);
    }
    
    
  
    fprintf(stdout,  "\n\n\n Random image (200x100x150 3D 16bit), 1 MPI proc, 1 thread, lambda \n" );  fflush(stdout);
    system("mpirun -np 1 ./disccofan -g 1,1,1 --inprefix check_files/Check16b --intype fits -o test --threads 1 -v off");
    if( access("out.fits", F_OK ) != -1){
      //     system("./polyImage compare check_files/Check16b_check fits out fits 1");
      compare_files("check_files/Check16b_check","fits", "out", "fits", 1, NULL, NULL);

      system("rm *fits");
    } else {
      fprintf(stdout, RED "Error, output files not written \n" RESET);  fflush(stdout);
      exit(0);
    }
  
    
    fprintf(stdout,  "\n Random image (200x100x150 3D 16bit), 32 MPI proc, 1 thread \n" );  fflush(stdout);
    system("mpirun  -np 32 ./disccofan -g 4,2,4 --inprefix check_files/Check16b --intype fits -o test --threads 1 -v off");
    if( access("out.fits", F_OK ) != -1){
      //  system("./polyImage compare check_files/Check16b_check fits out fits 1");
        compare_files("check_files/Check16b_check","fits", "out", "fits", 1, NULL, NULL);

      system("rm  *fits");
    } else {
      fprintf(stdout, RED "Error, output files not written \n" RESET);  fflush(stdout);
      exit(0);
    }
      
    fprintf(stdout,  "\n\n\n Random image (200x100x150 3D 16bit), 4 MPI proc, 3 threads \n" );  fflush(stdout);
    system("mpirun  -np 8 ./disccofan -g 2,2,2 --inprefix check_files/Check16b --intype fits -o test --threads 3 -v off");
    if( access("out.fits", F_OK ) != -1){
      // system("./polyImage compare check_files/Check16b_check fits out fits 1");
      compare_files("check_files/Check16b_check","fits", "out", "fits", 1, NULL, NULL);

      system("rm *fits");
    } else {
      fprintf(stdout, RED "Error, output files not written \n" RESET);  fflush(stdout);
      exit(0);
    }

    fprintf(stdout,  "\n\n\n Random image (200x100x150 3D 16bit), 4 MPI proc, 3 threads, connectivity 26 \n" );  fflush(stdout);
    system("mpirun  -np 8 ./disccofan -g 2,2,2 --inprefix check_files/Check16b --intype fits -o test --threads 3 -v off -c 26");
    if( access("out.fits", F_OK ) != -1){
      // system("./polyImage compare check_files/Check16b_check26 fits out fits 1");
      compare_files("check_files/Check16b_check26","fits", "out", "fits", 1, NULL, NULL);

      system("rm *fits");
    } else {
      fprintf(stdout, RED "Error, output files not written \n" RESET);  fflush(stdout);
      exit(0);
    }
  }

  if(choice == CSL || choice == ALL){
    fprintf(stdout, CYN "\n\n Checking CSL segmentation: start\n" RESET);  fflush(stdout);
    fprintf(stdout,  "\n\n\n Lena image (512x512 2D 8bit), 8 MPI proc, 4 thread, lvec.txt \n" );  fflush(stdout);
    system("mpirun  -np 8 ./disccofan -g 4,2,1 --inprefix check_files/lena512 --intype pgm --threads 4 -v off -o csl --outtype fits");
    if( access("out-C.fits", F_OK ) != -1 && access( "out-L.fits", F_OK ) != -1 && access( "out-S.fits", F_OK ) != -1) {
      //   system("./polyImage compare check_files/C pgm out-C fits 1");
      //   system("./polyImage compare check_files/L pgm out-L fits 1");
   //   system("./polyImage compare check_files/S pgm out-S fits 1");
      compare_files("check_files/C","pgm", "out-C", "fits", 1, NULL, NULL);
      compare_files("check_files/L","pgm", "out-L", "fits", 1, NULL, NULL);
      compare_files("check_files/S","pgm", "out-S", "fits", 1, NULL, NULL);

      system("rm *fits");
    } else {
      fprintf(stdout, RED "Error, output files not written \n" RESET);  fflush(stdout);
      exit(0);
    }

      
    fprintf(stdout,  "\n Random image (200x100x150 3D 16bit), 8 MPI proc, 3 threads, lvec.txt \n" );  fflush(stdout);
    system("mpirun  -np 8 ./disccofan -g 4,2,1 --inprefix check_files/Check16b --intype fits  --threads 3 -v off -o csl");
    if( access("out-C.fits", F_OK ) != -1 && access( "out-L.fits", F_OK ) != -1 && access( "out-S.fits", F_OK ) != -1) {
      //system("./polyImage compare check_files/Check-C fits out-C fits 1");
      //system("./polyImage compare check_files/Check-L fits out-L fits 1");
      // system("./polyImage compare check_files/Check-S fits out-S fits 1");
      compare_files("check_files/Check-C","fits", "out-C", "fits", 1, NULL, NULL);
      compare_files("check_files/Check-L","fits", "out-L", "fits", 1, NULL, NULL);
      compare_files("check_files/Check-S","fits", "out-S", "fits", 1, NULL, NULL);
      system("rm *fits");
    } else {
      fprintf(stdout, RED "Error, output files not written \n" RESET);  fflush(stdout);
      exit(0);
    }
  }

  if(choice == PATTERN || choice == ALL){
    fprintf(stdout, CYN "\n\n Checking Pattern spectrum: start\n" RESET);  fflush(stdout);
    fprintf(stdout,  "\n\n\n Lena image (512x512 2D 8bit), 8 MPI proc, 4 thread, lvec.txt \n" );  fflush(stdout);
    system("mpirun  -np 8 ./disccofan -g 4,2,1 --inprefix check_files/lena512 --intype pgm --threads 4 -v off -o pattern --lvec check_files/lvec.txt");
    if( access("pattern.txt", F_OK ) != -1) {
      // system("./polyImage compare check_files/lenapatt txt pattern txt 1");
      compare_files("check_files/lenapatt","txt", "pattern", "txt", 1, NULL, NULL);

      system("rm *txt");
    } else {
      fprintf(stdout, RED "Error, output files not written \n" RESET);  fflush(stdout);
      exit(0);
    }

      
    fprintf(stdout,  "\n Random image (200x100x150 3D 16bit), 8 MPI proc, 3 threads, lvec.txt \n" );  fflush(stdout);
    system("mpirun -np 8 ./disccofan -g 4,2,1 --inprefix check_files/Check16b --intype fits  --threads 3 -v off -o pattern --lvec check_files/lvec.txt");
    if( access("pattern.txt", F_OK ) != -1) {
      // system("./polyImage compare check_files/Checkpat txt pattern txt 1");
      compare_files("check_files/Checkpat","txt", "pattern", "txt", 1, NULL, NULL);

      system("rm *txt");
    } else {
      fprintf(stdout, RED "Error, output files not written \n" RESET);  fflush(stdout);
      exit(0);
    }
  }

  fprintf(stdout, GRN "\n\n\nCheck done \n" RESET);  fflush(stdout);
  return 0;
  
}
