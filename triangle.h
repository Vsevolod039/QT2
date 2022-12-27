# include <vector>
# include <math.h>
# include <stdio.h>
# include <stdlib.h>

class triangle_3D
{
  public:
    std::vector<double> *A = nullptr;
    std::vector<double> *B = nullptr;
    std::vector<double> *C = nullptr;
    std::vector<double> M = {0.0, 0.0, 0.0};
    double key = 0.0;
  public:
    friend class general;
    void set_triangle (std::vector<double> *sA, std::vector<double> *sB, std::vector<double> *sC);
    void mass_center ();
    void key_calculate (double x, double y, double z);
    int operator< (const triangle_3D &t) const;
    void print ();

    std::vector<double> getA() const;
};

class Point
{
  public:
    std::vector<double> p = {0.0, 0.0, 0.0};
  public:
    void set_point (double a, double b, double c)
    {
      p [0] = a, p [1] = b, p [2] = c;
    }

};
