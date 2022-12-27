//# include <sys/sysinfo.h>
# include <sys/types.h>
# include <sched.h>
# include <unistd.h>
# include <sys/resource.h>
# include <sys/time.h>
# include <pthread.h>

enum class task_from_GUI
{
  NO_TASK = 0,
  FUNCTION_CHANGED = 1,
  N_CHANGED = 2,
};

enum class status
{
  FREE = 0,
  BUSY = 1,
  READY = 2,
};

class general_data
{
  public :
    status s = status::FREE; // status of calculative threads
    task_from_GUI q = task_from_GUI::NO_TASK;
    int nx = 0;
    int ny = 0;
    int iter = 0;
    bool first = true;
    double left_upper_x = 0.0;
    double left_upper_y = 0.0;
    double right_down_x = 0.0;
    double right_down_y = 0.0;
    double (*f) (double, double);
    double splash = 0.0;
    double x_s = 0.0;
    double y_s = 0.0;
    // for GUI and CALC
    double *x;
    // for GUI
    //double *answer;
    //for CALC
    double *b;
    double *u;
    double *v;
    double *r;
    double *A;
    double *buf;
    int *I;
    pthread_mutex_t all = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_t calc = PTHREAD_MUTEX_INITIALIZER; // ONLY FOR CALCULATION !!!
    pthread_cond_t c_in = PTHREAD_COND_INITIALIZER; // ONLY FOR CALCULATION !!!
  public:
    ~general_data ()
    {
      delete [] x;
      //delete [] answer;
      /*delete [] b;
      delete [] u;
      delete [] v;
      delete [] r;
      delete [] A;
      delete [] buf;
      delete [] I;*/
    }

};
