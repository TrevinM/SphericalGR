#ifndef SCALARPOTENTIAL_H
#define SCALARPOTENTIAL_H

#include <iostream>
#include <fstream>

using namespace std;

class QuadraticPotential
{
  public:
    QuadraticPotential()
    {
        input_file = "QuadraticPotential_Input";
        ifstream infile;
        infile.open(input_file);
        if (!input_file)
        {
            cerr << "Couldn't open input for potential :(" << endl;
        }
        char buf[100], c;
        infile.get(buf, 100, '='); infile.get(c); infile >> mass;
        infile.close();

        cout << "Creating quadratic potential with mass m = " << mass << endl;
    }

    inline double V(double sf)
    {
        return 0.5 * mass * mass * sf * sf;
    }
    inline double dVdsf(double sf)
    {
        return mass*mass*sf;
    }

  private:
    double mass;
    char* input_file;
};

#endif
