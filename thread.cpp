# include "thread.h"

void * thread_func (void * ptr)
{
  args * a = (args * )ptr;
  while (1)
    {
      pthread_mutex_lock (&(a->arg->calc));
      while (a->arg->q == task_from_GUI::NO_TASK)
        {
          pthread_cond_wait (&(a->arg->c_in), &(a->arg->calc));
        }
      pthread_mutex_unlock (&(a->arg->calc));

      if (a->arg->q == task_from_GUI::N_CHANGED || a->arg->q == task_from_GUI::FUNCTION_CHANGED )
        {
          pthread_mutex_lock (&(a->arg->all));
          a->arg->s = status::BUSY;
          pthread_mutex_unlock (&(a->arg->all));
          a->nx = a->arg->nx;
          a->ny = a->arg->ny;
          a->func = a->arg->f;
          a->splash = a->arg->splash;
          a->x_s = a->arg->x_s;
          a->y_s = a->arg->y_s;
          int nx = a->nx;
          int ny = a->ny;
          if (a->k == 0)
            {
              pthread_mutex_lock (&(a->arg->all));
              a->arg->first = false;
              a->arg->b = new double [(nx + 1) * (ny + 1)];
              a->arg->x = new double [(nx + 1) * (ny + 1)];
              a->arg->r = new double [(nx + 1) * (ny + 1)];
              a->arg->u = new double [(nx + 1) * (ny + 1)];
              a->arg->v = new double [(nx + 1) * (ny + 1)];
              a->arg->A = new double [(nx + 1) * (ny + 1) + 1 + a->get_len (nx, ny)];
              a->arg->I = new int [(nx + 1) * (ny + 1) + 1 + a->get_len (nx, ny)];
              a->arg->buf = new double [a->p];
              a->allocate_msr_matrix (a->arg->A, a->arg->b, a->arg->I, a->nx, a->ny);
              pthread_mutex_unlock (&(a->arg->all));
            }
          a->left_upper_x = a->arg->left_upper_x;
          a->left_upper_y = a->arg->left_upper_y;
          a->right_down_x = a->arg->right_down_x;
          a->right_down_y = a->arg->right_down_y;

          if (a->k == 0)
            {
              for (int u = 0; u < (nx + 1) * (ny + 1); u ++)
                {
                  a->arg->x[u] = 0.0;
                  a->arg->b[u] = 0.0;
                  a->arg->r[u] = 0.0;
                  a->arg->u[u] = 0.0;
                  a->arg->v[u] = 0.0;
                }
            }
          reduce_sum (a->p, 0, 1);
          a->build_msr_matrix (a->arg->A, a->arg->I, a->nx, a->ny, a->p, a->k);
          a->build_colomn (a->arg->b, a->nx, a->ny, a->p, a->k);
          int res = a->minimal_residual_msr_matrix_full (a->arg->A, a->arg->I, (a->nx + 1) * (a->ny + 1), a->arg->b,
                                                           a->arg->x, a->arg->r, a->arg->u,
                                                           a->arg->v, a->eps, 1000, a->p, a->k, a->arg->buf);
          if (res < 0)
            {
              printf ("  No solution\n");
            }
          reduce_sum (a->p, 0, 1);
          pthread_mutex_lock (&(a->arg->all));
          a->arg->q = task_from_GUI::NO_TASK;
          a->arg->s = status::READY;
          a->arg->iter = a->iter;
          pthread_mutex_unlock (&(a->arg->all));
        }
    }
  return 0;
}
