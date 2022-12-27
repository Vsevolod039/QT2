# ifndef _GNU_SOURCE
# define _GNU_SOURCE
# endif
# include <stdio.h>
# include <stdlib.h>
//# include <sys/sysinfo.h>
# include <sys/types.h>
# include <sched.h>
# include <unistd.h>
# include <sys/resource.h>
# include <sys/time.h>
# include <pthread.h>
# include <math.h>
# include "general.h"
# include "func.h"

void reduce_sum (int, int *, int);

class MSR_matrix
{
  public:
    int N = 0;
    int LEN = 0; // Len of not diagonal part;
    double (*func) (double, double);

    double left_upper_x = 0.0;
    double left_upper_y = 0.0;
    double right_down_x = 0.0;
    double right_down_y = 0.0;

    int iter = 0;

    double splash = 0.0;
    double x_s = 0.0;
    double y_s = 0.0;
  public:
    int get_len (int nx, int ny);
    void ij2l (int nx, int ny, int i, int j, int &l);
    void l2ij (int nx, int ny, int &i, int &j, int l);

    int get_off_diag_num (int nx, int ny, int l);
    int get_off_diag_index (int nx, int ny, int l, int *J, double *a_diag, double *a);

    int allocate_msr_matrix (double *A, double *b, int *I, int nx, int ny);
    void build_msr_matrix (double *A, int *I, int nx, int ny, int p, int k);

    void mult_msr_matrix_vector (double *A, int *I, int n, double *x, double *y, int p, int k);
    void apply_preconditioner_msr_matrix (double *A, int n, double *y, double *x, int p, int k);
    double scalar_prod (int n, double *x, double *y, int p, int k, double * buf);
    void mult_add_vector (int n, double *x, double *y, double t, int p, int k);

    int minimal_residual_msr_matrix (double *A, int *I, int n, double *b, double *x, double *r, double *u, double * v,
                                                 double eps, int max_it, int p, int k, double *buf);
    int minimal_residual_msr_matrix_full (double *A, int *I, int n, double *b, double *x, double *r, double *u,
                                                      double *v, double eps, int max_it, int p, int k, double * buf);

    double integrate (int nx, int ny, int l);
    void build_colomn (double *b, int nx, int ny, int p, int k);

    void get_begin_end (int p, int n, int k, int &begin, int &end);

    void set_func (int id);
    void print_matrix (double *A, int *I, int nx, int ny);
    void print_colomn (double *b);
    void normalization (double *A, double *b, int nx, int ny, int p, int k);
    double f (double x0, double y0);

};


class args: public MSR_matrix
{
  public:
    int p = 0;
    int k = 0;
    int nx = 0;
    int ny = 0;
    double eps = 0.0;
    general_data *arg = nullptr;
  public:
    void set_arg (int sp, int sk, int snx, int sny, double seps, general_data *sarg)
    {
      p = sp;
      k = sk;
      nx = snx;
      ny = sny;
      eps = seps;
      arg = sarg;
    }
};
