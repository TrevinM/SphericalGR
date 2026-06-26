#include <iostream>
#include <functional>

class ClassName
{
public:
  double add(double a, double b);
  typedef double (ClassName::*Combiner)(double, double);
  double intermediate(double a, double b, Combiner);
  double combiner(double a, double b);
};

double ClassName::add(double a, double b)
{
  return a+b;
}

double ClassName::intermediate(double a, double b, Combiner func)
{
  return (this->*func)(a, b);
}

double ClassName::combiner(double a, double b)
{
  return intermediate(a, b, &ClassName::add);
}

int main()
{
  ClassName OBJ;
  std::cout << OBJ.combiner(12, 10) << std::endl;
}

