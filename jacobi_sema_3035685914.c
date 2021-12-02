/* 
// filename: jacobi_sema_3035685914.c
// student name: Mak Tsz Shing
// student number: 3035685914
// development platform: vi editor via workbench2
// remark: 
*/
#define _GNU_SOURCE

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>              /* For timing */
#include <sys/time.h>            /* For timing */
#include <sys/resource.h>
#include <string.h>
#include <pthread.h>
#include <stdbool.h>
#include <semaphore.h>


/****************Global****************************/

#define MAX(a,b) ((a)>(b)?(a):(b))
#define EPSILON 0.001            /* Termination condition */

char *filename;                  /* File name of output file */

/* Grid size */
int M = 200;                     /* Number of rows */
int N = 200;                     /* Number of cols */
long max_its = 1000000;          /* Maximum iterations, a safe bound to avoid infinite loop */
double final_diff;               /* Temperature difference between iterations at the end */

/* Thread count */
int thr_count = 2;

/* shared variables between threads */
/*************************************************************/
double** u;                   /* Previous temperatures */
double** w;                   /* New temperatures */


// (1) Add your variables here

typedef struct room room;
struct room {
	int Mi;
	int Mj;
	int N;
};
sem_t diff_mutex, sig_wrt, sig_read;
int diff_i = 0;
bool sig = false;
int sig_readcnt = 0;
bool finish;
/**************************************************************/

int main (int argc, char *argv[])
{
   int      its;                 /* Iterations to converge */
   double   elapsed;             /* Execution time */
   struct timeval stime, etime;  /* Start and end times */
   struct rusage usage;

   void allocate_2d_array (int, int, double ***);
   void initialize_array (double ***);
   void print_solution (char *, double **);
   int  find_steady_state (void);

   /* For convenience of other problem size testing */
   if ((argc == 1) || (argc == 4)) {
      if (argc == 4) {
         M = atoi(argv[1]);
         N = atoi(argv[2]);
         thr_count = atoi(argv[3]);
      } // Otherwise use default grid and thread size
   } else {
     printf("Usage: %s [ <rows> <cols> <threads ]>\n", argv[0]);
     exit(-1);
   }

   printf("Problem size: M=%d, N=%d\nThread count: T=%d\n", M, N, thr_count);

   /* Create the output file */
   filename = argv[0];
   sprintf(filename, "%s.dat", filename);

   allocate_2d_array (M, N, &u);
   allocate_2d_array (M, N, &w);
   initialize_array (&u);
   initialize_array (&w);

   gettimeofday (&stime, NULL);
   its = find_steady_state();
   gettimeofday (&etime, NULL);

   elapsed = ((etime.tv_sec*1000000+etime.tv_usec)-(stime.tv_sec*1000000+stime.tv_usec))/1000000.0;

   printf("Converged after %d iterations with error: %8.6f.\n", its, final_diff);
   printf("Elapsed time = %8.4f sec.\n", elapsed);

   getrusage(RUSAGE_SELF, &usage);
   printf("Program completed - user: %.4f s, system: %.4f s\n",
      (usage.ru_utime.tv_sec + usage.ru_utime.tv_usec/1000000.0),
    (usage.ru_stime.tv_sec + usage.ru_stime.tv_usec/1000000.0));
   printf("no. of context switches: vol %ld, invol %ld\n\n",
  		  usage.ru_nvcsw, usage.ru_nivcsw);

   print_solution (filename, w);
}

/* Allocate two-dimensional array. */
void allocate_2d_array (int r, int c, double ***a)
{
   double *storage;
   int     i;
   storage = (double *) malloc (r * c * sizeof(double));
   *a = (double **) malloc (r * sizeof(double *));
   for (i = 0; i < r; i++)
      (*a)[i] = &storage[i * c];
}

/* Set initial and boundary conditions */
void initialize_array (double ***u)
{
   int i, j;

   /* Set initial values and boundary conditions */
   for (i = 0; i < M; i++) {
      for (j = 0; j < N; j++)
         (*u)[i][j] = 25.0;      /* Room temperature */
      (*u)[i][0] = 0.0;
      (*u)[i][N-1] = 0.0;
   }

   for (j = 0; j < N; j++) {
      (*u)[0][j] = 0.0;
      (*u)[M-1][j] = 1000.0;     /* Heat source */
   }
}

/* Print solution to standard output or a file */
void print_solution (char *filename, double **u)
{
   int i, j;
   char sep;
   FILE *outfile;

   if (!filename) { /* if no filename specified, print on screen */
      sep = '\t';   /* tab added for easier view */
      outfile = stdout;
   } else {
      sep = '\n';   /* for gnuplot format */
      outfile = fopen(filename,"w");
      if (outfile == NULL) {
         printf("Can't open output file.");
         exit(-1);
      }
   }

   /* Print the solution array */
   for (i = 0; i < M; i++) {
      for (j = 0; j < N; j++)
         fprintf (outfile, "%6.2f%c", u[i][j], sep);
      fprintf(outfile, "\n"); /* Empty line for gnuplot */
   }
   if (outfile != stdout)
      fclose(outfile);

}

/* Entry function of the worker threads */
void *thr_func(void *arg) {

// (2) Add the worker's logic here
    struct rusage usage;

    //Identify which set of rows to compute
	room *rm = arg;
	int mi = rm->Mi;
	int mj = rm->Mj;
	int n = rm->N;

	double diff;
	int its;
    bool wake;

	for (its = 1; its < max_its; its++) {
		diff = 0.0;
        //Compute the temp of all points in its set
		for (int i = mi; i <= mj; i++) {
			for (int j = 1; j < n-1; j++) {
				(w)[i][j] = 0.25 * ((u)[i-1][j] + (u)[i+1][j] + (u)[i][j-1] + (u)[i][j+1]);
                if (fabs((w)[i][j] - (u)[i][j]) > diff)
                    diff = fabs((w)[i][j] - (u)[i][j]);
			}
        }

		//write to diff
        sem_wait(&diff_mutex);
        final_diff = MAX(final_diff, diff);
        diff_i++;
        sem_post(&diff_mutex);

        sem_wait(&sig_read);
        wake = sig;

		if (!wake) {
			break;
		}
	}
	
    //Get its running statistics and pass it to master
    getrusage(RUSAGE_SELF, &usage);
    struct rusage *rtn = malloc(sizeof(struct rusage));
    *rtn = usage;
    //terminate
    pthread_exit(rtn);
}


int find_steady_state (void)
{

// (3) Implement the thread creation and the main control logic here
    struct rusage master_usage, *thr_usage;

	sem_init(&diff_mutex, 0, 1);
	sem_init(&sig_wrt, 0, 1);
	sem_init(&sig_read, 0, 0);

    double **temp;
    int temp_count;
    void *retval;
    
    double diff = 0.0;
    bool signal = false;
	int thr_m[thr_count];
	int sum = 0;
	int div = (M-2) / thr_count;
	int rmd = (M-2) % thr_count;
	for (int i = 0; i < thr_count; i++) {
		if (rmd > 0) {
			thr_m[i] = div + 1;
			rmd--;
		}
		else {
			thr_m[i] = div;
		}
	}

    //Create N workers; each with a set of rows
	pthread_t thread[thr_count];
	for (int i = 0; i < thr_count; i++) {
		room *rm = malloc(sizeof(room));
		rm->Mi = sum + 1;
		rm->Mj = thr_m[i] + sum;
		sum += thr_m[i];
		rm->N = N;
		int s = pthread_create(&thread[i], NULL, thr_func, (void*) rm);
		if (s != 0) {
			perror("pthread_create: ");
			exit(-1);
		}
	}
	
	while (!finish) {
        if (diff_i == thr_count) {
            diff_i = 0;

            sem_wait(&diff_mutex);
            //printf("final_diff = %f\n", final_diff);
            temp = u;
            u = w;
            w = temp;
            temp_count++;
            diff = final_diff;
            final_diff = 0.0;
            sem_post(&diff_mutex);

            //check final_diff against EPSILON
            if (diff <= EPSILON) {
                signal = false;
                finish = true;
            }
            else {
                signal = true;
            }

            sem_wait(&sig_wrt);
            sig = signal;
            sem_post(&sig_wrt);

            for (int i = 0; i < thr_count; i++) {
                sem_post(&sig_read);
            }
        }
    }
    for (int j = 0; j < thr_count; j++) {
        pthread_join(thread[j], &retval);
        thr_usage = (struct rusage*)retval;
        printf("Thread %d has completed - user: %.4f s, system: %.4f s\n", j,
        (thr_usage->ru_utime.tv_sec + thr_usage->ru_utime.tv_usec/1000000.0),
        (thr_usage->ru_stime.tv_sec + thr_usage->ru_stime.tv_usec/1000000.0));
    }
    
    getrusage(RUSAGE_SELF, &master_usage);
    printf("find_stedy_state - user: %.4f s, system: %.4f s\n",
        (master_usage.ru_utime.tv_sec + master_usage.ru_utime.tv_usec/1000000.0),
        (master_usage.ru_stime.tv_sec + master_usage.ru_stime.tv_usec/1000000.0));

    return temp_count;
}
