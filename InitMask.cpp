// created by Gabriel Balestre, Stéphanie Lyon, Nathan Paillou, Théo Robushi
// supervised by Annabelle Collin and Heloise Beaugendre
// 2017/18

#include "InitMask.h"
#include "Util.h"
#include "LevelSet_v.h"
#include <iostream>
#include <fstream>

using namespace std;

InitMask::InitMask(){};

std::vector<std::vector<double>> InitMask::BuildMaskAndRedistancing(int rows, int cols, std::string nameMaskRedist)
{
  std::vector<std::vector<double>> M_v;
  std::ifstream file(nameMaskRedist.data());
  if (file)
  {
    M_v=readVTKFile(M_v, nameMaskRedist);

    if ( (M_v.size() != rows) || (M_v[0].size() != cols) )
    {
      std::cout << "The initialization " << nameMaskRedist << " does not correspond to the image (size problem)." << std::endl;
      abort();
    }
  }
  else
  {
    bool is_ok = true;
    if (rows > cols)
    {
      is_ok = false;
      int temp = cols;
      cols = rows;
      rows = temp;
    }

    M_v.resize(rows);
    for (size_t i = 0; i < rows; i++)
    {
      M_v[i].resize(cols);
    }

    int num_of_circles = 2;
    const double radius = rows/(3.*num_of_circles+1.);
    int num_cols = floor(cols/(3.*radius));
    double dist = (cols-num_cols*2.*radius)/(1.0+num_cols);

    // #pragma acc parallel loop reduction(+:M_v)
    std:: cout << M_v.size() << " " << M_v[0].size()<< std::endl;

    for(int i=0; i<rows; i++)
    {
      for(int j=0; j<cols; j++)
      {
        M_v[i][j]=0;
        for (int k=0; k<num_cols; k++)
        {
          for (int l=0; l<num_of_circles; l++)
          {
            double cy((dist+radius)+k*radius+2*k*dist);
            double cx((2+3*l)*radius);
            if  (sqrt((i-cx)*(i-cx)+(j-cy)*(j-cy)) < radius)
            {
              M_v[i][j] = 1;
            }
            else
            {
              M_v[i][j]=-1;
            }
          }
        }
      }
    }

    std::cout << "Redistanciation en cours ... " << std::endl;

    LevelSet_v lv(M_v);
    lv.redistancing_v(200);
    saveVTKFile(M_v, nameMaskRedist);

  }
  return M_v;
}
