
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
sem_t diff_empty, diff_full, diff_mutex, sig_empty, sig_full, sig_mutex;
bool sig = false;
int thr_finish_count;
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

	for (its = 1; its < max_its; its++) {
		diff = 0.0;
        //Compute the temp of all points in its set
		for (int i = mi; i <= mj; i++) {
			for (int j = 1; j < n-1; j++) {
				w[i][j] = 0.25 * (u[i-1][j] + u[i+1][j] + u[i][j-1] + u[i][j+1]);
                //Find the max temp diff in this set
				diff = MAX(fabs(w[i][j] - u[i][j]), diff);
			}
        }

		//Signal master that to test the diff
		sem_wait(&diff_empty);
		sem_wait(&diff_mutex);
		final_diff = diff;
		sem_post(&diff_mutex);
		sem_post(&diff_full);

		//consume wake_sig
		bool wake;
		sem_wait(&sig_full);
		sem_wait(&sig_mutex);
		wake = sig;
		sig = false;
		sem_post(&sig_mutex);
		sem_post(&sig_empty);

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

	sem_init(&diff_empty, 0, 1);
	sem_init(&diff_full, 0, 0);
	sem_init(&diff_mutex, 0, 1);
	sem_init(&sig_empty, 0, 1);
	sem_init(&sig_full, 0, 0);
	sem_init(&sig_mutex, 0, 1);

    double **temp;
    int temp_count;
    void *retval;
    double diff = 0.0;
	int thr_m[thr_count];
	int sum = 0;
	int div = (M-1) / thr_count;
	int rmd = (M-1) % thr_count;
	for (int i = 0; i < thr_count; i++) {
		if (rmd == 0) {
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
	
	while (thr_finish_count < thr_count) {
		//consume final_diff and produce wake_sig
		sem_wait(&diff_full);
		sem_wait(&diff_mutex);
		sem_wait(&sig_empty);
		sem_wait(&sig_mutex);

        //swap u, w
        temp = u;
        u = w;
        w = temp;

        temp_count++;

        //check final_diff against EPSILON
		if (final_diff <= EPSILON) {
			sig = false;
			thr_finish_count++;
		}
		else {
			sig = true;
		}
		sem_post(&sig_mutex);
		sem_post(&sig_full);
        diff = final_diff;
		final_diff = 0.0;
		sem_post(&diff_mutex);
		sem_post(&diff_empty);
	}
	
    /*
	void *retval;
	int sum_its = 0;
	for (int j = 0; j < thr_count; j++) {
		pthread_join(thread[j], &retval);
		sum_its += (int)retval;
	}
    */

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

    final_diff = diff;

    return temp_count;
}
