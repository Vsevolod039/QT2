# include "func.h"

double functions::f0 (double x, double y)
{
  (void)x, (void)y;
  return 1.0;
}

double functions::f1 (double x, double y)
{
  (void)y;
  return x;
}

double functions::f2 (double x, double y)
{
  (void)x;
  return y;
}

double functions::f3 (double x, double y)
{
  return x + y;
}

double functions::f4 (double x, double y)
{
  return sqrt (x * x + y * y);
}

double functions::f5 (double x, double y)
{
  return x * x + y * y;
}

double functions::f6 (double x, double y)
{
  return exp (x * x - y * y);
}

double functions::f7 (double x, double y)
{
  return 1 / (25 * (x * x + y * y) + 1);
}
