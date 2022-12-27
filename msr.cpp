# include "msr.h"

double MSR_matrix::f (double x0, double y0)
{
  if (fabs (x0 - x_s) < 1e-6 && fabs (y0 - y_s) < 1e-6)
    return func (x0, y0) + splash;
  return func (x0, y0);
}

int MSR_matrix::get_len (int nx, int ny)
{
  LEN = 6 * (nx - 1) * (ny - 1) + 4 * (2 * (nx - 1) + 2 * (ny - 1)) + 3 * 2 + 2 * 2;
  return LEN;
}

void MSR_matrix::ij2l (int nx, int ny, int i, int j, int &l)
{
  (void)ny;
  l = i + j * (nx + 1);
}

void MSR_matrix::l2ij (int nx, int ny, int &i, int &j, int l)
{
  (void)ny;
  j = l / (nx + 1);
  i = l - j * (nx + 1);
}

int MSR_matrix::get_off_diag_index (int nx, int ny, int l, int *J, double *a_diag, double *a)
{
  int i = 0, j = 0;
  l2ij (nx, ny, i, j, l);
  double hx = (right_down_x - left_upper_x) / nx, hy = (left_upper_y - right_down_y) / ny;
  if (i > 0 && i < nx && j > 0 && j < ny)
    {
      *a_diag = (6.0 / 12.0) * hx * hy;
      // 0
      // neighbors (i + 1, j), (i, j - 1), (i - 1, j + 1), (i - 1, j), (i, j + 1), (i + 1, j + 1)
      ij2l (nx, ny, i + 1, j + 0, J [0]);
      ij2l (nx, ny, i + 0, j - 1, J [1]);
      ij2l (nx, ny, i - 1, j - 1, J [2]);
      ij2l (nx, ny, i - 1, j + 0, J [3]);
      ij2l (nx, ny, i + 0, j + 1, J [4]);
      ij2l (nx, ny, i + 1, j + 1, J [5]);
      for (int w = 0; w < 6; w ++)
        a [w] = (1.0 / 12.0) * hx * hy;
      return 6;
    }
  if (i == 0 && j > 0 && j < ny)
    {
      // neighbors (i + 1, j), (i, j - 1), (i, j + 1), (i + 1, j + 1)
      // 1
      *a_diag = (3.0 / 12.0) * hx * hy;
      ij2l (nx, ny, i + 1, j + 0, J [0]);
      ij2l (nx, ny, i + 0, j - 1, J [1]);
      ij2l (nx, ny, i + 0, j + 1, J [2]);
      ij2l (nx, ny, i + 1, j + 1, J [3]);
      a [0] = (1.0 / 12.0) * hx * hy;
      a [1] = (1.0 / 24.0) * hx * hy;
      a [2] = (1.0 / 24.0) * hx * hy;
      a [3] = (1.0 / 12.0) * hx * hy;
      return 4;
    }
  if (i == nx && j > 0 && j < ny)
    {
      // neighbors (i, j - 1), (i - 1, j - 1), (i - 1, j), (i, j + 1)
      // 2
      *a_diag = (3.0 / 12.0) * hx * hy;
      ij2l (nx, ny, i + 0, j - 1, J [0]);
      ij2l (nx, ny, i - 1, j - 1, J [1]);
      ij2l (nx, ny, i - 1, j + 0, J [2]);
      ij2l (nx, ny, i + 0, j + 1, J [3]);
      a [0] = (1.0 / 24.0) * hx * hy;
      a [1] = (1.0 / 12.0) * hx * hy;
      a [2] = (1.0 / 12.0) * hx * hy;
      a [3] = (1.0 / 24.0) * hx * hy;
      return 4;
    }
  if (i > 0 && i < nx && j == 0)
    {
      // neighbors (i + 1, j), (i - 1, j), (i, j + 1), (i + 1, j + 1)
      // 3
      *a_diag = (3.0 / 12.0) * hx * hy;
      ij2l (nx, ny, i + 1, j + 0, J [0]);
      ij2l (nx, ny, i - 1, j + 0, J [1]);
      ij2l (nx, ny, i + 0, j + 1, J [2]);
      ij2l (nx, ny, i + 1, j + 1, J [3]);
      a [0] = (1.0 / 24.0) * hx * hy;
      a [1] = (1.0 / 24.0) * hx * hy;
      a [2] = (1.0 / 12.0) * hx * hy;
      a [3] = (1.0 / 12.0) * hx * hy;
      return 4;
    }
  if (i > 0 && i < nx && j == ny)
    {
      // neighbors (i + 1, j), (i, j - 1), (i - 1, j - 1), (i - 1, j)
      // 4
      *a_diag = (3.0 / 12.0) * hx * hy;
      ij2l (nx, ny, i + 1, j + 0, J [0]);
      ij2l (nx, ny, i + 0, j - 1, J [1]);
      ij2l (nx, ny, i - 1, j - 1, J [2]);
      ij2l (nx, ny, i - 1, j + 0, J [3]);
      a [0] = (1.0 / 24.0) * hx * hy;
      a [1] = (1.0 / 12.0) * hx * hy;
      a [2] = (1.0 / 12.0) * hx * hy;
      a [3] = (1.0 / 24.0) * hx * hy;
      return 4;
    }
  if (i == 0 && j == 0)
    {
      // neighbors (i + 1, j), (i, j + 1), (i + 1, j + 1)
      // 5
      *a_diag = (2.0 / 12.0) * hx * hy;
      ij2l (nx, ny, i + 1, j + 0, J [0]);
      ij2l (nx, ny, i + 0, j + 1, J [1]);
      ij2l (nx, ny, i + 1, j + 1, J [2]);
      a [0] = (1.0 / 24.0) * hx * hy;
      a [1] = (1.0 / 24.0) * hx * hy;
      a [2] = (1.0 / 12.0) * hx * hy;
      return 3;
    }
  if (i == nx && j == ny)
    {
      // neighbors (i, j - 1), (i - 1, j - 1), (i - 1, j)
      // 6
      *a_diag = (2.0 / 12.0) * hx * hy;
      ij2l (nx, ny, i + 0, j - 1, J [0]);
      ij2l (nx, ny, i - 1, j - 1, J [1]);
      ij2l (nx, ny, i - 1, j + 0, J [2]);
      a [0] = (1.0 / 24.0) * hx * hy;
      a [1] = (1.0 / 12.0) * hx * hy;
      a [2] = (1.0 / 24.0) * hx * hy;
      return 3;
    }
  if (i == 0 && j == ny)
    {
      // neighbors (i + 1, j), (i, j - 1)
      // 7
      *a_diag = (1.0 / 12.0) * hx * hy;
      ij2l (nx, ny, i + 1, j + 0, J [0]);
      ij2l (nx, ny, i + 0, j - 1, J [1]);
      a [0] = (1.0 / 24.0) * hx * hy;
      a [1] = (1.0 / 24.0) * hx * hy;
      return 2;
    }
  if (i == nx && j == 0)
    {
      // neighbors (i - 1, j), (i, j + 1)
      // 8
      *a_diag = (1.0 / 12.0) * hx * hy;
      ij2l (nx, ny, i - 1, j + 0, J [0]);
      ij2l (nx, ny, i + 0, j + 1, J [1]);
      a [0] = (1.0 / 24.0) * hx * hy;
      a [1] = (1.0 / 24.0) * hx * hy;
      return 2;
    }
  return -1;
}

int MSR_matrix::allocate_msr_matrix (double *A, double *b, int *I, int nx, int ny)
{
  (void)b;
  int i = 0, l = 0, s = 0;
  N = (nx + 1) * (ny + 1); // diag
  int nz = get_len (nx, ny);
  int len = N + 1 + nz;

  if (!A)
    return -1;
  if (!I)
    {
      delete [] A;
      return -2;
    }
  I [0] = N + 1;
  for (i = 1; i <= N; i ++)
    {
      l = get_off_diag_num (nx, ny, i - 1);
      I [i] = l + I [i - 1];
      s += l;
    }
  if (s != nz)
    {
      printf ("  ******************************\n");
      printf ("  **** ABORT (string 179),   ***\n");
      printf ("  *******************************\n");
      return -3;
    }
  if (I [N] != len)
    {
      printf ("  ******************************\n");
      printf ("  **** ABORT (string 186),   ***\n");
      printf ("  *******************************\n");
      return -4;
    }
  return 0;
}

void MSR_matrix::build_msr_matrix (double *A, int *I, int nx, int ny, int p, int k)
{
  int N = (nx + 1) * (ny + 1), sum = 0;
  int begin = 0, end = 0, i = 0;
  get_begin_end (p, N, k, begin, end);

  for (i = begin; i < end; i ++)
    {
     sum += get_off_diag_index (nx, ny, i, I + I [i], A + i, A + I [i]);
    }

  A [N] = 0.0;

  reduce_sum (p, &sum, 1);
  if (sum + N + 1 != I [N])
    {
      printf ("  ******************************\n");
      printf ("  **** ABORT (string 209), %d ***\n", k);
      printf ("  *******************************\n");
    }
}

int MSR_matrix::get_off_diag_num (int nx, int ny, int l)
{
  int i = 0, j = 0;
  l2ij (nx, ny, i, j, l);
  if (i > 0 && i < nx && j > 0 && j < ny) return 6;
  if ((i == 0 || i == nx) && (j > 0 && j < ny)) return 4;
  if ((i > 0 && i < nx) && (j == 0 || j == ny)) return 4;
  if ((i == 0 && j == 0) || (i == nx && j == ny)) return 3;
  if ((i == 0 && j == ny) || (i == nx && j == 0)) return 2;
  return -1;
}

void MSR_matrix::mult_msr_matrix_vector (double *A, int *I, int n, double *x, double *y, int p, int k)
{
  int i = 0, j = 0, J = 0;
  int begin = 0, end = 0, l = 0;
  double s = 0.0;
  get_begin_end (p, n, k, begin, end);

  for (i = begin; i < end; i ++)
    {
      s = A [i] * x [i];
      l = I [i + 1] - I [i];
      J = I [i];
      for (j = 0; j < l; j ++)
        s += A [J + j] * x [I [J + j]];
      y [i] = s;
    }
  reduce_sum (p, 0, 1);
}

void MSR_matrix::apply_preconditioner_msr_matrix (double *A, int n, double *y, double *x, int p, int k)
{
  int begin = 0, end = 0, i = 0;
  get_begin_end (p, n, k, begin, end);

  for (i = begin; i < end; i ++)
    x [i] = y [i] / A [i];
  reduce_sum (p, 0, 1);
}

double MSR_matrix::scalar_prod (int n, double *x, double *y, int p, int k, double * buf)
{
  int begin = 0, end = 0, i = 0;
  double s = 0.0;
  //begin = k * n, begin /= p;
  //end = (k + 1) * n, end /= p;
  get_begin_end (p, n, k, begin, end);

  for (i = begin; i < end; i ++)
    s += x [i] * y [i];

  buf [k] = s;
  reduce_sum (p, 0, 1);
  s = 0.0;

  for (i = 0; i < p; i ++)
    s += buf [i];

  reduce_sum (p, 0, 1);
  return s;
}

void MSR_matrix::mult_add_vector (int n, double *x, double *y, double t, int p, int k)
{
  int begin = 0, end = 0, i = 0;
  get_begin_end (p, n, k, begin, end);

  for (i = begin; i < end; i ++)
    x [i] -= t * y [i];
  reduce_sum (p, 0, 1);
}

int MSR_matrix::minimal_residual_msr_matrix (double *A, int *I, int n, double *b, double *x, double *r, double *u, double * v,
                                             double eps, int max_it, int p, int k, double *buf)
{
  int it = 0;
  double C1 = 0.0, C2 = 0.0, tau = 0.0;
  double b_norm2 = scalar_prod (n, b, b, p, k, buf);
  double prec = eps * eps * b_norm2;
  mult_msr_matrix_vector (A, I, n, x, r, p, k);
  mult_add_vector (n, r, b, 1.0, p, k);

  for (it = 0; it < max_it; it ++)
    {
      apply_preconditioner_msr_matrix(A, n, r, v, p, k);
      mult_msr_matrix_vector(A, I, n, v, u, p, k);
      C1 = scalar_prod (n, u, r, p, k, buf);
      C2 = scalar_prod (n, u, u, p, k, buf);
      if (C1 < prec || C2 < prec)
        break;
      tau = C1 / C2;
      mult_add_vector (n, x, v, tau, p, k);
      mult_add_vector (n, r, u, tau, p, k);
      reduce_sum (p, 0, 1);
    }

  reduce_sum (p, 0, 1);
  //printf ("from %d = it = %d\n", k, it);
  if (it >= max_it)
    return -1;

  reduce_sum (p, 0, 1);
  return it;
}

int MSR_matrix::minimal_residual_msr_matrix_full (double *A, int *I, int n, double *b, double *x, double *r, double *u,
                                                  double *v, double eps, int max_it, int p, int k, double * buf)
{
  int step = 50, ret = 0, it = 0;
  for (it = 0; it < max_it; it += step)
    {
      ret += minimal_residual_msr_matrix (A, I, n, b, x, r, u, v, eps, step, p, k, buf);
      reduce_sum (p, 0, 1);
      if (ret >= 0)
        break;
    }
  iter = ret;
  reduce_sum (p, 0, 1);
  if (it >= max_it)
    return -1;
  return 0;
}

void MSR_matrix::get_begin_end (int p, int n, int k, int &begin, int &end)
{
  if (n % p == 0)
    {
      begin = k * (n / p);
      end = (k + 1) * (n / p);
    }
  else
    {
      if (k < n % p)
        {
          begin = k * (n / p) + k;
          end = (k + 1) * (n / p) + k + 1;
        }
      else
        {
          begin = n % p + k * (n / p);
          end = n % p + (k + 1)* (n / p);
        }
    }
}

double MSR_matrix::integrate (int nx, int ny, int l)
{
  int i = 0, j = 0;
  l2ij(nx, ny, i, j, l);
  double hx = (right_down_x - left_upper_x) / nx, hy = (left_upper_y - right_down_y) / ny;
  double x0 = left_upper_x + i * hx, y0 = right_down_y + j * hy;

  double l1 = 0.0, l2 = 0.0, l3 = 0.0, l4 = 0.0;

  if (i > 0 && i < nx && j > 0 && j < ny)
    {
      l1 = f (x0, y0);
      l2 = f (x0 + hx / 2, y0) + f (x0 + hx / 2, y0 + hy / 2) + f (x0, y0 + hy / 2) +
           f (x0 - hx / 2, y0) + f (x0 - hx / 2, y0 - hy / 2) + f (x0, y0 - hy / 2);
      l3 = f (x0 + hx, y0 + hy / 2) +  f (x0 + hx / 2, y0 + hy) + f (x0 - hx / 2, y0 + hy / 2) +
           f (x0 - hx, y0 - hy / 2) + f (x0 - hx / 2, y0 - hy) + f (x0 + hx / 2, y0 - hy / 2);
      l4 = f (x0 + hx, y0) + f (x0 + hx, y0 + hy) + f (x0, y0 + hy) +
           f (x0 - hx, y0) + f (x0 - hx, y0 - hy) + f (x0, y0 - hy);

      //printf ("1  l1 = %lf, l2 = %lf, l3 = %lf, l4 = %lf\n", l1, l2, l3, l4);
      return  (36 * l1 + 20 * l2 + 4 * l3 + 2 * l4) / 192;
    }
  if (i == 0 && j > 0 && j < ny)
    {
      l1 = f (x0, y0);
      l2 = 2 * f (x0 + hx / 2, y0) + 2 * f (x0 + hx / 2, y0 + hy / 2) + f (x0, y0 + hy / 2) +
           f (x0, y0 - hy / 2);
      l3 = f (x0 + hx, y0 + hy / 2) + f (x0 + hx / 2, y0 + hy) + f (x0 + hx / 2, y0 - hy / 2);
      l4 = 2 * f (x0 + hx, y0) + 2 * f (x0 + hx, y0 + hy) + f (x0, y0 + hy) + f (x0, y0 - hy);

      //printf ("2  l1 = %lf, l2 = %lf, l3 = %lf, l4 = %lf\n", l1, l2, l3, l4);
      return  (18 * l1 + 10 * l2 + 4 * l3 + 1 * l4) / 192;
    }
  if (i == nx && j > 0 && j < ny)
    {
      l1 = f (x0, y0);
      l2 = f (x0, y0 + hy / 2) + 2 * f (x0 - hx / 2, y0) + 2 * f (x0 - hx / 2, y0 - hy / 2) +
           f (x0, y0 - hy / 2);
      l3 = f (x0 - hx / 2, y0 + hy / 2) + f (x0 - hx, y0 - hy / 2) + f (x0 - hx / 2, y0 - hy);
      l4 = f (x0, y0 + hy) + 2 * f (x0 - hx, y0) + 2 * f (x0 - hx, y0 - hy) + f (x0, y0 - hy);

      //printf ("3  l1 = %lf, l2 = %lf, l3 = %lf, l4 = %lf\n", l1, l2, l3, l4);
      return  (18 * l1 + 10 * l2 + 4 * l3 + 1 * l4) / 192;
    }
  if (i > 0 && i < nx && j == 0)
    {
      l1 = f (x0, y0);
      l2 = f (x0 + hx / 2, y0) + 2 * f (x0 + hx / 2, y0 + hy / 2) + 2 * f (x0, y0 + hy / 2) +
           f (x0 - hx / 2, y0);
      l3 = f (x0 + hx, y0 + hy / 2) + f (x0 + hx / 2, y0 + hy) + f (x0 - hx / 2, y0 + hy / 2);
      l4 = f (x0 + hx, y0) + 2 * f (x0 + hx, y0 + hy) + 2 * f (x0, y0 + hy) + f (x0 - hx, y0);

      //printf ("4  l1 = %lf, l2 = %lf, l3 = %lf, l4 = %lf\n", l1, l2, l3, l4);
      return  (18 * l1 + 10 * l2 + 4 * l3 + 1 * l4) / 192;
    }
  if (i > 0 && i < nx && j == ny)
    {
      l1 = f (x0, y0);
      l2 = f (x0 + hx / 2, y0) + f (x0 - hx / 2, y0) + 2 * f (x0 - hx / 2, y0 - hy / 2) +
           2 * f (x0, y0 - hy / 2);
      l3 = f (x0 - hx, y0 - hy / 2) + f (x0 - hx / 2, y0 - hy) + f (x0 + hx / 2, y0 - hy / 2);
      l4 = f (x0 - hx, y0) + 2 * f (x0 - hx, y0 - hy) + 2 * f (x0, y0 - hy) + f (x0 + hx, y0);

      //printf ("5  l1 = %lf, l2 = %lf, l3 = %lf, l4 = %lf\n", l1, l2, l3, l4);
      return  (18 * l1 + 10 * l2 + 4 * l3 + 1 * l4) / 192;
    }
  if (i == 0 && j == 0)
    {
      l1 = f (x0, y0);
      l2 = f (x0 + hx / 2, y0) + 2 * f (x0 + hx / 2, y0 + hy / 2) + f (x0, y0 + hy / 2);
      l3 = f (x0 + hx, y0 + hy / 2) + f (x0 + hx / 2, y0 + hy);
      l4 = f (x0 + hx, y0) + 2 * f (x0 + hx, y0 + hy) + f (x0, y0 + hy);

      //printf ("6  l1 = %lf, l2 = %lf, l3 = %lf, l4 = %lf\n", l1, l2, l3, l4);
      return  (12 * l1 + 10 * l2 + 4 * l3 + 1 * l4) / 192;
    }
  if (i == nx && j == ny)
    {
      l1 = f (x0, y0);
      l2 = f (x0 - hx / 2, y0) + 2 * f (x0 - hx / 2, y0 - hy / 2) + f (x0, y0 - hy / 2);
      l3 = f (x0 - hx, y0 - hy / 2) + f (x0 - hx / 2, y0 - hy); //+ f (x0 + hx / 2, y0 - hy / 2);
      l4 = f (x0 - hx, y0) + 2 * f (x0 - hx, y0 - hy) + f (x0, y0 - hy);

      //printf ("7  l1 = %lf, l2 = %lf, l3 = %lf, l4 = %lf\n", l1, l2, l3, l4);
      return  (12 * l1 + 10 * l2 + 4 * l3 + 1 * l4) / 192;
    }
  if (i == 0 && j == ny)
    {
      l1 = f (x0, y0);
      l2 = f (x0 + hx / 2, y0) + f (x0, y0 - hy / 2);
      l3 = f (x0 + hx / 2, y0 - hy / 2);
      l4 = f (x0 + hx, y0) + f (x0, y0 - hy);

      //printf ("8  l1 = %lf, l2 = %lf, l3 = %lf, l4 = %lf\n", l1, l2, l3, l4);
      return  (6 * l1 + 10 * l2 + 4 * l3 + 1 * l4) / 192;
    }
  if (i == nx && j == 0)
    {
      l1 = f (x0, y0);
      l2 = f (x0, y0 + hy / 2) + f (x0 - hx / 2, y0);
      l3 = f (x0 - hx / 2, y0 + hy / 2);
      l4 = f (x0, y0 + hy) + f (x0 - hx, y0);

      //printf ("9  l1 = %lf, l2 = %lf, l3 = %lf, l4 = %lf\n", l1, l2, l3, l4);
      return  (6 * l1 + 10 * l2 + 4 * l3 + 1 * l4) / 192;
    }
  return -1;
}

void MSR_matrix::build_colomn (double *b, int nx, int ny, int p, int k)
{
  int begin = 0, end = 0, N = (nx + 1) * (ny + 1), i = 0;
  get_begin_end(p, N, k, begin, end);
  double hx = (right_down_x - left_upper_x) / nx, hy = (left_upper_y - right_down_y) / ny;
  for (i = begin; i < end; i ++)
    {
      b [i] = (integrate (nx, ny, i)) * hx * hy;
    }
  reduce_sum (p, 0, 1);
}

void MSR_matrix::print_matrix (double *A, int *I, int nx, int ny)
{
  int i = 0, k = 0;
  printf ("  ************  My  MSR  matrix  ************  \n");
  printf ("  size = %d\n", N);
  printf ("  Pattern : LEN = %d, N = %d\n", LEN, N);
  for (i = 0; i < N; i ++)
    {
      printf (" %.3lf ", A [i]);
      k ++;
    }
  printf ("NO\n");
  k ++;
  for (i = 0; i < N; i ++)
    {
      int l = get_off_diag_num (nx, ny, i);
      printf ("Scalar product for %d ======= \n", i);
      for (int j = 0; j < l; j ++)
        {
          printf (" %d ", k);
          printf (" %.3lf (%d)",  A[k], I [k]);
          k ++;
        }
      printf ("\n");
    }
  printf ("\n");
}

void MSR_matrix::print_colomn (double *b)
{
  for (int i = 0; i < N; i ++)
    {
      printf ("%10.3lf\n", b [i]);
    }
}

void MSR_matrix::normalization (double *A, double *b, int nx, int ny, int p, int k)
{
  int i = 0;
  (void)p, (void)k;
  double hx = (right_down_x - left_upper_x) / nx, hy = (left_upper_y - right_down_y) / ny;
  for (i = 0; i < N + LEN + 1; i ++)
    A [i] = A [i] * hx * hy;
  for (i = 0; i < N; i ++)
    b [i] = b [i] * hx * hy;

}

void reduce_sum (int p, int * a, int n)
{
  static pthread_mutex_t m = PTHREAD_MUTEX_INITIALIZER;
  static pthread_cond_t c_in = PTHREAD_COND_INITIALIZER;
  static pthread_cond_t c_out = PTHREAD_COND_INITIALIZER;
  static int t_in = 0;
  static int t_out = 0;
  static int * p_a = 0;
  int i = 0;
  if (p <= 1)
    return;
  pthread_mutex_lock (&m);
  if (!p_a)
    {
      p_a = a;
    }
  else
    {
      for (i = 0; i < n; i ++)
        p_a [i] += a[i];
    }
  t_in ++;
  if (t_in >= p)
    {
      t_out = 0;
      pthread_cond_broadcast (&c_in);
    }
  else
    {
      while (t_in < p)
        {
          pthread_cond_wait (& c_in, &m);
        }
    }
  if (p_a != a)
    {
      for (i = 0; i < n; i ++)
        a [i] = p_a [i];
    }
  t_out ++;
  if (t_out >= p)
    {
      t_in = 0;
      p_a = 0;
      pthread_cond_broadcast (&c_out);
    }
  else
    {
      while (t_out < p)
        {
          pthread_cond_wait (&c_out, &m);
        }
    }
  pthread_mutex_unlock (&m);
}

