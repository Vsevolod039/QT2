#include <QPainter>
#include <stdio.h>

#include "general.h"

# define LEN 1234

# define EPS 1e-14

double X0 = 1;
double Y0 = 1;
double Z0 = 1;

General::General (QWidget *parent)
  : QWidget (parent)
{
  func_id = 0;
}

void General::set_arg (general_data *a)
{
  arg = a;
}

double General::func (double x0, double y0)
{
  if (fabs (x0 - x_s) < 1e-6 && fabs (y0 - y_s) < 1e-6)
    return f (x0, y0) + splash;
  return f (x0, y0);
}

void General::timer ()
{
  status s_loc = status::FREE;
  pthread_mutex_lock (&(arg->all));
  s_loc = arg->s;
  pthread_mutex_unlock (&(arg->all));

  if (s_loc == status::READY)
    {
      delete [] x;
      delete [] arg->A;
      delete [] arg->I;
      delete [] arg->buf;
      delete [] arg->b;
      delete [] arg->r;
      delete [] arg->u;
      delete [] arg->v;

      //printf ("  nx = %d, ny = %d\n", nx, ny);
      //if (v == visual::FUNCTION)  printf ("  F_max = %e,  F_min = %e\n", F_max, F_min);
      //if (v == visual::APPROXIMATION)  printf ("  AP_max = %e,  AP_min = %e\n", AP_max, AP_min);
      //if (v == visual::RESIDUAL)  printf ("  R_max = %e,  R_min = %e\n", R_max, R_min);
      //printf ("  Iterations = %d\n", arg->iter);
      x = arg->x;
      m = moving::WAS_MOVED;
      nx = nx_new;
      ny = ny_new;
      pthread_mutex_lock (&(arg->all));
      arg->s = status::FREE;
      pthread_mutex_unlock (&(arg->all));
      update ();
    }
}

int General::parse_command_line (int argc, char **argv)
{
  fname = argv [1];
  FILE * fp = fopen (fname, "r");
  if (argc != 7)
    {
      printf ("  User : ");
      printf ("not enought arguments!\n");
      printf ("  %s  fname.txt  nx  ny  k  eps  p !\n", argv [0]);
      fclose (fp);
      return -1;
    }

  if (   !fp
      || sscanf (argv [2], "%d", &nx) != 1
      || sscanf (argv [3], "%d", &ny) != 1
      || sscanf (argv [4], "%d", &k) != 1
      || sscanf (argv [5], "%lf", &eps) != 1
      || sscanf (argv [6], "%d", &p) != 1)
    {
      printf ("  User : ");
      printf ("%s  fname.txt  nx  ny  k  eps  p !\n", argv [0]);
      fclose (fp);
      return -2;
    }

  arg->nx = nx;
  arg->ny = ny;

  nx_new = nx;
  ny_new = ny;

  int res = read_file (fp);
  if (res)
    {
      printf ("  User : ");
      printf ("not correct input file!\n");
      fclose (fp);
      return -3;
    }

  arg->left_upper_x = left_upper_x;
  arg->left_upper_y = left_upper_y;
  arg->right_down_x = right_down_x;
  arg->right_down_y = right_down_y;

  func_id = (k - 1) % 8;
  change_func ();

  fclose (fp);
  return 0;
}

int General::read_file (FILE * fp)
{
  char buf [LEN] = {0};
  int param = 0;

  while (fgets (buf, LEN, fp))
    {
      int i = 0;
      while (buf [i] == ' ')
        i ++;

      if (buf [i] == '#')
        continue;

      if (param == 0)
        {
          if (sscanf (buf + i, "%lf", &left_upper_x) == 1)
            {
              while (buf [i] != ' ')
                i ++;

              if (sscanf (buf + i + 1, "%lf", &left_upper_y) == 1)
                param = 1;
            }
        }
      else
        {
          if (sscanf (buf + i, "%lf", &right_down_x) == 1)
            {
              while (buf [i] != ' ')
                i ++;

              if (sscanf (buf + i + 1, "%lf", &right_down_y) == 1)
                return 0;
            }
        }
    }
  return -1;
}

void General::print_argument ()
{
  printf ("  INT argument : nx = %d  ny = %d  k = %d  p = %d\n", nx, ny, k, p);
  printf ("  DOUBLE argument : eps = %lf| Г (%.3lf, %.3lf) _| (%.3lf, %.3lf)\n",
          eps, left_upper_x, left_upper_y, right_down_x, right_down_y);
  printf ("  STRING argument : %s\n", fname);
}

void General::go_through_area ()
{
  int i = 0, j = 0;
  double dx = (right_down_x - left_upper_x) / Nx;
  double dy = (left_upper_y - right_down_y) / Ny;

  first = false;
  if (v == visual::FUNCTION)
    {
      for (i = 0; i < (Ny + 1); i ++)
        {
          for (j = 0; j < Nx + 1; j ++)
            {
              double x = left_upper_x + dx * j;
              double y = left_upper_y - dy * i;
              P [j + i * (Ny + 1)] = {x, y, func (x, y)};
            }
        }
    }

  if (v == visual::APPROXIMATION)
    {
      for (i = 0; i < (Ny + 1); i ++)
        {
          for (j = 0; j < Nx + 1; j ++)
            {
              double x = left_upper_x + dx * j;
              double y = left_upper_y - dy * i;
              P [j + i * (Ny + 1)] = {x, y, app (x, y)};
            }
        }
    }

  if (v == visual::RESIDUAL)
    {
      for (i = 0; i < (Ny + 1); i ++)
        {
          for (j = 0; j < Nx + 1; j ++)
            {
              double x = left_upper_x + dx * j;
              double y = left_upper_y - dy * i;
              P [j + i * (Ny + 1)] = {x, y, fabs (func (x, y) - app (x, y))};
            }
        }
    }

  for (i = 0; i < 2 * Ny; i += 2)
    {
      for (j = 0; j < Nx; j ++)
        {
          triangle_array [Ny * (i + 0) + j].set_triangle (P + j + i / 2 * (Ny + 1),
                                                          P + j + 1 + i / 2 * (Ny + 1),
                                                          P + j + (i / 2 + 1) * (Ny + 1));
          triangle_array [Ny * (i + 0) + j].mass_center ();
          triangle_array [Ny * (i + 0) + j].key_calculate (X0, Y0, Z0);

          triangle_array [Ny * (i + 1) + j].set_triangle (P + j + 1 + i / 2 * (Ny + 1),
                                                          P + j + (i / 2 + 1) * (Ny + 1),
                                                          P + j + 1 + (i / 2 + 1) * (Ny + 1));
          triangle_array [Ny * (i + 1) + j].mass_center ();
          triangle_array [Ny * (i + 1) + j].key_calculate (X0, Y0, Z0);
        }
    }
}
// z1 = 1, z2 = z3 = 0, calculate in x0, y0
double General::basic (double x1, double y1, double x2, double y2, double x3, double y3, double x0, double y0)
{
  return ((x0 - x1) * (y3 - y2) +
          (y0 - y1) * (x2 - x3)) / ((x3 - x1) * (y2 - y1) - (x2 -x1) * (y3 - y1)) + 1;

}

int binar_search_for_x (double x0, double a, double b, double n)
{
  if (x0 <= a) return 0;
  if (x0 >= b) return n;

  int i = 0, j = n, z = 0;
  double h = (b - a) / n;

  while (i != j)
    {
      z = (i + j) / 2;
      if (x0 <= a + z * h)
        {
          j = z;
        }
      else
        {
          i = z + 1;
        }
    }
  return i;
}

int binar_search_for_y (double x0, double a, double b, double n)
{
  if (x0 <= a) return 0;
  if (x0 >= b) return n;

  int i = 0, j = n, z = 0;
  double h = (b - a) / n;

  while (i != j)
    {
      z = (i + j) / 2;
      if (x0 <= a + z * h)
        {
          j = z;
        }
      else
        {
          i = z + 1;
        }
    }
  return i - 1;
}


double General::app (double x0, double y0)
{
  int i = 0, j = 0;
  double dx = (right_down_x - left_upper_x) / nx;
  double dy = (left_upper_y - right_down_y) / ny;

  i = binar_search_for_x (x0, left_upper_x, right_down_x, nx);
  j = ny - binar_search_for_y (y0, right_down_y, left_upper_y, ny);

  if (fabs (y0 - (left_upper_y - dy * j)) < EPS && fabs (x0 - (left_upper_x + dx * i)) < EPS)
    {
      return x [(ny - j) * (nx + 1) + i];
    }

  if (fabs (y0 - left_upper_y) < EPS)
    {
      double x1 = left_upper_x + dx * (i - 1);
      double y1 = left_upper_y - dy * (j + 0);

      double x2 = left_upper_x + dx * (i + 0);
      double y2 = left_upper_y - dy * (j + 0);

      double x3 = left_upper_x + dx * (i - 1);
      double y3 = left_upper_y - dy * (j + 1);

      double ans = x [(ny - (j + 0)) * (nx + 1) + i + 0] * basic (x2, y2, x1, y1, x3, x3, x0, y0) +
                   x [(ny - (j + 0)) * (nx + 1) + i - 1] * basic (x1, y1, x2, y2, x3, y3, x0, y0);
      return ans;
    }

  if (fabs (y0 - right_down_y) < EPS)
    {
      double x1 = left_upper_x + dx * (i - 1);
      double y1 = left_upper_y - dy * (j + 0);

      double x2 = left_upper_x + dx * (i + 0);
      double y2 = left_upper_y - dy * (j + 0);

      double x3 = left_upper_x + dx * (i + 0);
      double y3 = left_upper_y - dy * (j - 1);

      double ans = x [(ny - (j + 0)) * (nx + 1) + i + 0] * basic (x2, y2, x1, y1, x3, x3, x0, y0) +
                   x [(ny - (j + 0)) * (nx + 1) + i - 1] * basic (x1, y1, x2, y2, x3, y3, x0, y0);
      return ans;
    }

  if (fabs (x0 - left_upper_x) < EPS)
    {
      double x1 = left_upper_x + dx * (i + 0);
      double y1 = left_upper_y - dy * (j - 1);

      double x2 = left_upper_x + dx * (i + 0);
      double y2 = left_upper_y - dy * (j + 0);

      double x3 = left_upper_x + dx * (i + 1);
      double y3 = left_upper_y - dy * (j - 1);

      double ans = x [(ny - (j + 0)) * (nx + 1) + i + 0] * basic (x2, y2, x1, y1, x3, x3, x0, y0) +
                   x [(ny - (j - 1)) * (nx + 1) + i + 0] * basic (x1, y1, x2, y2, x3, y3, x0, y0);
      return ans;
    }

  if (fabs (x0 - right_down_x) < EPS)
    {
      double x1 = left_upper_x + dx * (i + 0);
      double y1 = left_upper_y - dy * (j - 1);

      double x2 = left_upper_x + dx * (i + 0);
      double y2 = left_upper_y - dy * (j + 0);

      double x3 = left_upper_x + dx * (i - 1);
      double y3 = left_upper_y - dy * (j + 0);

      double ans = x [(ny - (j + 0)) * (nx + 1) + i + 0] * basic (x2, y2, x1, y1, x3, x3, x0, y0) +
                   x [(ny - (j - 1)) * (nx + 1) + i + 0] * basic (x1, y1, x2, y2, x3, y3, x0, y0);
      return ans;
    }

  double x1 = left_upper_x + dx * (i - 1);
  double y1 = left_upper_y - dy * (j - 1);

  double x2 = left_upper_x + dx * (i + 0);
  double y2 = left_upper_y - dy * (j - 1);

  double x3 = left_upper_x + dx * (i - 1);
  double y3 = left_upper_y - dy * (j + 0);

  double x4 = left_upper_x + dx * (i + 0);
  double y4 = left_upper_y - dy * (j + 0);

  double ro_1 = (x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0);
  double ro_2 = (x4 - x0) * (x4 - x0) + (y4 - y0) * (y4 - y0);

  if (ro_1 < ro_2)
    {
      double ans = x [(ny - (j - 1)) * (nx + 1) + i - 1] * basic (x1, y1, x2, y2, x3, y3, x0, y0) +
             x [(ny - (j - 1)) * (nx + 1) + i + 0] * basic (x2, y2, x1, y1, x3, y3, x0, y0) +
             x [(ny - (j + 0)) * (nx + 1) + i - 1] * basic (x3, y3, x1, y1, x2, y2, x0, y0);
      return ans;
    }
  else
    {
      double ans = x [(ny - (j + 0)) * (nx + 1) + i + 0] * basic (x4, y4, x2, y2, x3, y3, x0, y0) +
             x [(ny - (j - 1)) * (nx + 1) + i + 0] * basic (x2, y2, x4, y4, x3, y3, x0, y0) +
             x [(ny - (j + 0)) * (nx + 1) + i - 1] * basic (x3, y3, x4, y4, x2, y2, x0, y0);
      return ans;
    }

  printf ("  ******************************\n");
  printf ("  **** ABORT (string 306),   ***\n");
  printf ("  *******************************\n");
}

int comp (const void *val1, const void *val2)
{
  return (*((const triangle_3D *)val1) < *((const triangle_3D *)val2));
}

void General::sort_triangle ()
{
  qsort (triangle_array, 2 * Nx * Ny, sizeof (triangle_3D), comp);
}


void General::set_func()
{
  f = functions::f5;
}

void General::find_real_max_min ()
{
  // For function
  F_max = -1e10, F_min = 1e10;
  double value = 0.0;
  triangle_3D *T;
  for (int i = 0; i < 2 * Nx * Ny; i ++)
    {
      T = triangle_array + i;
      value = func ((*T->A) [0], (*T->A) [1]);
      if (F_max < value)
        F_max = value;

      if (F_min > value)
        F_min = value;
    }
  // For approximation
  AP_max = -1e10, AP_min = 1e10;
  value = 0.0;

  for (int i = 0; i < 2 * Nx * Ny; i ++)
    {
      T = triangle_array + i;
      value = app ((*T->A) [0], (*T->A) [1]);

      if (AP_max < value)
        AP_max = value;

      if (AP_min > value)
        AP_min = value;
    }

  R_max = -1e10, R_min = 1e10;
  value = 0.0;

  for (int i = 0; i < 2 * Nx * Ny; i ++)
    {
      T = triangle_array + i;
      value = (*T->A) [2];//fabs (func ((T->M) [0], (T->M) [1]) - app ((T->M) [0], (T->M) [1]));

      if (R_max < value)
        R_max = value;

      if (R_min > value)
        R_min = value;

      value = (*T->B) [2];
      if (R_max < value)
        R_max = value;

      if (R_min > value)
        R_min = value;

      value = (*T->C) [2];

      if (R_max < value)
        R_max = value;

      if (R_min > value)
        R_min = value;

    }
}

void General::find_min_max ()
{
  max = {0.0, 0.0, -1e10};
  min = {0.0, 0.0, 1e10};
  triangle_3D *T;
  for (int i = 0; i < 2 * Nx * Ny; i ++)
    {
      T = triangle_array + i;

      if (get_V (max [0], max [1], max [2]) < get_V ((*T->A) [0], (*T->A) [1], (*T->A) [2]))
        max = *(triangle_array [i]).A;

      if (get_V (max [0], max [1], max [2]) < get_V ((*T->B) [0], (*T->B) [1], (*T->B) [2]))
        max = *(triangle_array [i]).B;

      if (get_V (max [0], max [1], max [2]) < get_V ((*T->C) [0], (*T->C) [1], (*T->C) [2]))
        max = *(triangle_array [i]).C;

      if (get_V (min [0], min [1], min [2]) > get_V ((*T->A) [0], (*T->A) [1], (*T->A) [2]))
        min = *(triangle_array [i]).A;

      if (get_V (min [0], min [1], min [2]) > get_V ((*T->B) [0], (*T->B) [1], (*T->B) [2]))
        min = *(triangle_array [i]).B;

      if (get_V (min [0], min [1], min [2]) > get_V ((*T->C) [0], (*T->C) [1], (*T->C) [2]))
        min = *(triangle_array [i]).C;
    }
}

void General::increase_nx_and_ny ()
{
  status s_loc = status::FREE;
  pthread_mutex_lock (&(arg->all));
  s_loc = arg->s;
  pthread_mutex_unlock (&(arg->all));
  if (s_loc == status::BUSY)
    {
      printf ("  User : now threads is calculate\n");
      return;
    }
  nx_new *= 2;
  ny_new *= 2;
  m = moving::NOT_MOVE;
  pthread_mutex_lock (&(arg->all));
  arg->nx = nx_new;
  arg->ny = ny_new;
  arg->q = task_from_GUI::N_CHANGED;
  pthread_cond_broadcast (&(arg->c_in));
  pthread_mutex_unlock (&(arg->all));
}

void General::reduce_nx_and_ny ()
{
  status s_loc = status::FREE;
  pthread_mutex_lock (&(arg->all));
  s_loc = arg->s;
  pthread_mutex_unlock (&(arg->all));
  if (s_loc == status::BUSY)
    {
      printf ("  User : now threads is calculate\n");
      return;
    }
  if (nx / 2 > 3 && ny / 2 > 3)
    {
      nx_new /= 2;
      ny_new /= 2;
    }
  m = moving::NOT_MOVE;
  pthread_mutex_lock (&(arg->all));
  arg->nx = nx_new;
  arg->ny = ny_new;
  arg->q = task_from_GUI::N_CHANGED;
  pthread_cond_broadcast (&(arg->c_in));
  pthread_mutex_unlock (&(arg->all));
}

void General::zoom ()
{
  scale -= 50;
  update ();
}

void General::zoom_out ()
{
  scale += 50;
  update ();
}

void General::left_rotate ()
{
  fi += 3.1415 / 12;
  update ();
}

void General::right_rotate ()
{
  fi -= 3.1415 / 12;
  update ();
}

void General::increase_splash ()
{
  status s_loc = status::FREE;
  pthread_mutex_lock (&(arg->all));
  s_loc = arg->s;
  pthread_mutex_unlock (&(arg->all));
  if (s_loc == status::BUSY)
    {
      printf ("  User : now threads is calculate\n");
      return;
    }
  splash += 0.1 * F_max;
  m = moving::NOT_MOVE;
  x_s = (right_down_x - left_upper_x) / 2;
  y_s = (left_upper_y - right_down_y) / 2;
  pthread_mutex_lock (&(arg->all));
  arg->splash = splash;
  arg->x_s = (right_down_x - left_upper_x) / 2;
  arg->y_s = (left_upper_y - right_down_y) / 2;
  arg->q = task_from_GUI::FUNCTION_CHANGED;
  pthread_cond_broadcast (&(arg->c_in));
  pthread_mutex_unlock (&(arg->all));
}

void General::reduce_splash ()
{
  status s_loc = status::FREE;
  pthread_mutex_lock (&(arg->all));
  s_loc = arg->s;
  pthread_mutex_unlock (&(arg->all));
  if (s_loc == status::BUSY)
    {
      printf ("  User : now threads is calculate\n");
      return;
    }
  splash -= 0.1 * F_max;
  m = moving::NOT_MOVE;
  x_s = (right_down_x - left_upper_x) / 2;
  y_s = (left_upper_y - right_down_y) / 2;
  pthread_mutex_lock (&(arg->all));
  arg->splash = splash;
  arg->x_s = (right_down_x - left_upper_x) / 2;
  arg->y_s = (left_upper_y - right_down_y) / 2;
  arg->q = task_from_GUI::FUNCTION_CHANGED;
  pthread_cond_broadcast (&(arg->c_in));
  pthread_mutex_unlock (&(arg->all));
}


double General::get_U (double x, double y, double z)
{

  (void)z;
  return lambda * ((sin (fi) * (x + 0.5 * scale)  + cos (fi) * (y + 0.5 * scale))
                   - (cos (fi) * (x + 0.5 * scale) - sin (fi) * (y + 0.5 * scale)) / 2);
}

double General::get_V (double x, double y, double z)
{
  return lambda * ((z + 0.5 * scale) -
                   (sqrt (3) / 2) * (cos (fi) * (x + 0.5 * scale) - sin (fi) * (y + 0.5 * scale)));
}

void General::paintEvent(QPaintEvent *event)
{
  (void)event;
  QPainter painter (this);

  go_through_area ();
  sort_triangle ();
  find_min_max();
  if (fabs (splash - 0) < 1e-14 && (m == moving::WAS_MOVED)) find_real_max_min ();

  draw_function(&painter);
  int len = height () / 20;
  int sp_len = height () / 30;
  QString s_0 = QString ("[%1 ,  %2] x [%3 , %4]").arg(left_upper_x).arg(right_down_x).
      arg(right_down_y).arg(left_upper_y);
  painter.drawText (20, len, s_0);
  painter.drawText (20, len + sp_len, func_name);
  if (v == visual::FUNCTION)
    {
      QString s_1 = QString ("F_max = %1").arg(F_max);
      painter.drawText (20, len + 2 * sp_len, s_1);
      QString s_2 = QString ("F_min = %1").arg(F_min);
      painter.drawText (20, len + 3 * sp_len, s_2);
    }
  if (v == visual::APPROXIMATION)
    {
      QString s_1 = QString ("AP_max = %1").arg(AP_max);
      painter.drawText (20, len + 2 * sp_len, s_1);
      QString s_2 = QString ("AP_min = %1").arg(AP_min);
      painter.drawText (20, len + 3 * sp_len, s_2);
    }
  if (v == visual::RESIDUAL)
    {
      QString s_1 = QString ("R_max = %1").arg(R_max);
      painter.drawText (20, len + 2 * sp_len, s_1);
      QString s_2 = QString ("R_min = %1").arg(R_min);
      painter.drawText (20, len + 3 * sp_len, s_2);
    }
  QString s_3 = QString ("iterations = %1").arg(arg->iter);
  painter.drawText (20, len + 4 * sp_len, s_3);
  QString s_4 = QString ("nx = %1, ny = %2").arg(nx).arg(ny);
  painter.drawText (20, len + 5 * sp_len, s_4);
  printf ("  nx = %d, ny = %d, ", nx, ny);
  if (v == visual::FUNCTION)  printf ("  F_max = %.3e,  F_min = %.3e,", F_max, F_min);
  if (v == visual::APPROXIMATION)  printf ("  AP_max = %.3e,  AP_min = %.3e, ", AP_max, AP_min);
  if (v == visual::RESIDUAL)  printf ("  R_max = %.3e,  R_min = %.3e, ", R_max, R_min);
  printf ("  It = %d\n", arg->iter);
}

void General::draw_function(QPainter *painter)
{
  W_max = get_U (left_upper_x, left_upper_y, 0);
  W_min = get_U (right_down_x, right_down_y, 0);
  H_min = get_V (min [0], min [1], min [2]);
  H_max = get_V (max [0], max [1], max [2]);

  if (v == visual::FUNCTION)  painter->setBrush(QColor(255,165,0));
  if (v == visual::APPROXIMATION)  painter->setBrush(Qt::red);
  if (v == visual::RESIDUAL)  painter->setBrush(Qt::blue);

  for (int i = 0; i < 2 * Nx * Ny; i ++)
    {
      if (i >= 0)
        {
          triangle_3D *T = triangle_array + i;
          QPointF points [3] =
          {
            QPointF (U (get_U ((*T->A) [0], (*T->A) [1], (*T->A) [2])), V (get_V ((*T->A) [0], (*T->A) [1], (*T->A) [2]))),
            QPointF (U (get_U ((*T->B) [0], (*T->B) [1], (*T->B) [2])), V (get_V ((*T->B) [0], (*T->B) [1], (*T->B) [2]))),
            QPointF (U (get_U ((*T->C) [0], (*T->C) [1], (*T->C) [2])), V (get_V ((*T->C) [0], (*T->C) [1], (*T->C) [2]))),
          };
          painter->drawPolygon (points, 3);
        }
    }


}

double General::X (int u)
{
  return (u - 1) * (W_max - W_min) / (width () - scale) + W_min;
}

int General::V (double y)
{
  if (fabs (H_max - H_min) < 1e-6)
    {
      return height () / 2;
    }
  else
    {
      return height () - (int) (1 + (height () - scale) * (y - H_min) / (H_max - H_min));
    }
  return 0;
}

double General::Y (int v)
{
  return (v - 1) * (H_max - H_min) / (height () - scale) + H_min;
}

int General::U (double x)
{
  if (fabs (W_max - W_min) < 1e-6)
    {
      return 1;
    }
  else
    {
      return (int) (1 + (width () - scale) * (x - W_min) / (W_max - W_min));
    }
  return 0;
}

void General::change_func ()
{
  status s_loc = status::FREE;
  pthread_mutex_lock (&(arg->all));
  s_loc = arg->s;
  pthread_mutex_unlock (&(arg->all));
  if (s_loc == status::BUSY)
    {
      printf ("  User : now threads is calculate\n");
      return;
    }
  func_id = (func_id + 1) % 8;

  switch (func_id)
    {
      case 0:
        f = functions::f0;
        func_name = "function : f (x, y) = 1; func_id = 0";
        printf ("  =============== f (x, y) = 1 =======  func_id = 0\n");
        splash = 0.0;
        break;
      case 1:
        f = functions::f1;
        func_name = "function : f (x, y) = x; func_id = 1";
        printf ("  =============== f (x, y) = x =======  func_id = 1\n");
        splash = 0.0;
        break;
      case 2:
        f = functions::f2;
        func_name = "function : f (x, y) = y; func_id = 2";
        printf ("  =============== f (x, y) = y =======  func_id = 2\n");
        splash = 0.0;
        break;
      case 3:
        f = functions::f3;
        func_name = "function : f (x, y) = x + y; func_id = 3";
        printf ("  =============== f (x, y) = x + y =======  func_id = 3\n");
        splash = 0.0;
        break;
      case 4:
        f = functions::f4;
        func_name = "function : f (x, y) = sqrt (x^2 + y^2); func_id = 4";
        printf ("  =============== f (x, y) = sqrt (x^2 + y^2) =======  func_id = 4\n");
        splash = 0.0;
        break;
      case 5:
        f = functions::f5;
        func_name = "function : f (x, y) = x^2 + y^2; func_id = 5";
        printf ("  =============== f (x, y) = x^2 + y^2 =======  func_id = 5\n");
        splash = 0.0;
        break;
      case 6:
        f = functions::f6;
        func_name = "function : f (x, y) = exp (x^2 - y^2); func_id = 6";
        printf ("  =============== f (x, y) = exp (x^2 - y^2) =======  func_id = 6\n");
        splash = 0.0;
        break;
      case 7:
        f = functions::f7;
        func_name = "function : f (x, y) = 1 / (25 (x^2 + y^2) + 1); func_id = 7";
        printf ("  =============== f (x, y) = 1 / (25 (x^2 + y^2) + 1) =======  func_id = 7\n");
        splash = 0.0;
        break;
    }


  if (v == visual::FUNCTION)  printf ("  F_max = %e,  F_min = %e\n", F_max, F_min);
  if (v == visual::APPROXIMATION)  printf ("  AP_max = %e,  AP_min = %e\n", AP_max, AP_min);
  if (v == visual::RESIDUAL)  printf ("  R_max = %e,  R_min = %e\n", R_max, R_min);

  m = moving::NOT_MOVE;
  pthread_mutex_lock (&(arg->all));
  arg->f = f;
  arg->q = task_from_GUI::FUNCTION_CHANGED;
  arg->splash = 0.0;
  pthread_cond_broadcast (&(arg->c_in));
  pthread_mutex_unlock (&(arg->all));
}

void General::change_method ()
{
  status s_loc = status::FREE;
  pthread_mutex_lock (&(arg->all));
  s_loc = arg->s;
  pthread_mutex_unlock (&(arg->all));

  if (s_loc == status::BUSY)
    {
      printf ("  User : now threads is calculate\n");
      return;
    }
  switch (v)
    {
      case visual::FUNCTION:
        v = visual::APPROXIMATION;
        break;
      case visual::APPROXIMATION:
        v = visual::RESIDUAL;
        break;
      case visual::RESIDUAL:
        v = visual::FUNCTION;
      break;
    }
  /*if (v == visual::FUNCTION)  printf ("  F_max = %e,  F_min = %e\n", F_max, F_min);
  if (v == visual::APPROXIMATION)  printf ("  AP_max = %e,  AP_min = %e\n", AP_max, AP_min);
  if (v == visual::RESIDUAL)  printf ("  R_max = %e,  R_min = %e\n", R_max, R_min);*/

  update ();
}


