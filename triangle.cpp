# include "triangle.h"

std::vector<double> triangle_3D::getA() const
{
  return *A;
}

void triangle_3D::set_triangle (std::vector<double> *sA, std::vector<double> *sB, std::vector<double> *sC)
{
  A = sA, B = sB, C = sC;
}


void triangle_3D::mass_center()
{
  M [0] = ((*A) [0] + (*B) [0] + (*C) [0]) / 3.0;
  M [1] = ((*A) [1] + (*B) [1] + (*C) [1]) / 3.0;
  M [2] = ((*A) [2] + (*B) [2] + (*C) [2]) / 3.0;

}

void triangle_3D::key_calculate (double x, double y, double z)
{
  key = sqrt ((M [0] - x) * (M [0] - x) + (M [1] - y) * (M [1] - y) + (M [2] - z) * (M [2] - z));
}

int triangle_3D::operator< (const triangle_3D &t) const
{
  return (key < t.key);
}

void triangle_3D::print ()
{
  printf ("  >>>>>>>>>> A = (%.3lf,  %.3lf,  %.3lf)\n", (*A) [0], (*A) [1], (*A) [2]);
  printf ("  >>>>>>>>>> B = (%.3lf,  %.3lf,  %.3lf)\n", (*B) [0], (*B) [1], (*B) [2]);
  printf ("  >>>>>>>>>> C = (%.3lf,  %.3lf,  %.3lf)\n", (*C) [0], (*C) [1], (*C) [2]);
  printf ("  >>>>>>>>>> M = (%.3lf,  %.3lf,  %.3lf)\n", M [0], M [1], M [2]);
  printf ("  >>>>>>>>>> KEY = %lf\n", key);
}
