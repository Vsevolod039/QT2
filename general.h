#ifndef WINDOW_H
#define WINDOW_H
#include <vector>
# include "triangle.h"
# include "func.h"
# include "arg.h"
#include <QtWidgets/QtWidgets>
# define Nx 64
# define Ny 64

class triangle_3D;

enum class visual
{
  FUNCTION,
  APPROXIMATION,
  RESIDUAL,
};

enum class moving
{
  NOT_MOVE,
  WAS_MOVED,
};

class General : public QWidget
{
  Q_OBJECT

private:
  int nx = 0;
  int ny = 0;
  bool first = true;
  int nx_new = 0;
  int ny_new = 0;
  const char *fname = nullptr;
  int k = 0;
  double eps = 0.0;
  int p = 0;
  double fi = 0.0;
  double lambda = 1.0;
  int scale = 0;
  moving m = moving::NOT_MOVE;

  double left_upper_x = 0.0;
  double left_upper_y = 0.0;
  double right_down_x = 0.0;
  double right_down_y = 0.0;

  double (*f) (double, double);

  std::vector<double> max = {0.0, 0.0, -1e10};
  std::vector<double> min = {0.0, 0.0, 1e10};

  double W_max = 0.0;
  double H_max = 0.0;
  double W_min = 0.0;
  double H_min = 0.0;

  double F_max = 0.0;
  double F_min = 0.0;
  double AP_max = 0.0;
  double AP_min = 0.0;
  double R_max = 0.0;
  double R_min = 0.0;

  int func_id = 0;

  general_data *arg = nullptr;

  triangle_3D triangle_array [2 * Nx * Ny];
  std::vector<double> P [(Nx + 1) * (Ny + 1)];

  visual v = visual::FUNCTION;

  double * x = nullptr;
  const char *func_name = nullptr;

  double splash = 0.0;
  double x_s = 0.0;
  double y_s = 0.0;

public :
  General (QWidget *parent);
  int parse_command_line (int argc, char **argv);
  int read_file (FILE * fp);
  void print_argument ();
  double get_U (double x, double y, double z);
  double get_V (double x, double y, double z);
  void go_through_area ();
  void set_func ();
  void set_arg (general_data *);
  int get_p () {return p;}
  int get_nx () {return nx;}
  int get_ny () {return ny;}
  double get_eps () {return eps;}
  int get_left_upper_x () {return left_upper_x;}
  int get_left_upper_y () {return left_upper_y;}
  int get_right_down_x () {return right_down_x;}
  int get_right_down_y () {return right_down_y;}
  general_data* get_general_data () {return arg;}
  void set_general_data (general_data *sarg) { arg = sarg;}
  void sort_triangle ();
  void find_min_max ();
  void find_real_max_min ();
  void draw_function (QPainter *painter);
  double X (int u);
  double Y (int v);

  int U (double x);
  int V (double y);

  double app (double x0, double y0);
  // z1 = 1, z2 = z3 = 0, calculate in x0, y0
  double basic (double x1, double y1, double x2, double y2, double x3, double y3, double x0, double y0);
  double func (double x0, double y0);

public slots:
  void change_func ();
  void change_method ();
  void zoom ();
  void zoom_out ();
  void left_rotate ();
  void right_rotate ();
  void timer ();
  void increase_nx_and_ny ();
  void reduce_nx_and_ny ();
  void increase_splash ();
  void reduce_splash ();

protected:
  void paintEvent (QPaintEvent *event);
};


#endif
