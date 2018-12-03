// created by Gabriel Balestre, Stéphanie Lyon, Nathan Paillou, Théo Robushi
// supervised by Annabelle Collin and Heloise Beaugendre
// 2017/18

#include "InitMask.h"
#include "Util.h"

#include <iostream>
#include <fstream>

using namespace std;

InitMask::InitMask(){}

void InitMask::BuildMaskAndRedistancing(int rows, int cols, std::string nameMaskRedist)
{
  std::ifstream file(nameMaskRedist.data());
  if (file)
  {
    readVTKFile(_M, nameMaskRedist);
    if ( (_M.rows() != rows) || (_M.cols() != cols) )
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
    _M = (Eigen::MatrixXd::Ones(rows,cols)).array();

    _M *= (-1);
    int num_of_circles = 2;
    const double radius = rows/(3.*num_of_circles+1.);
    int num_cols = floor(cols/(3.*radius));
    double dist = (cols-num_cols*2.*radius)/(1.0+num_cols);
    for(int i=0; i<rows; i++)
    {
      for(int j=0; j<cols; j++)
      {
        for (int k=0; k<num_cols; k++)
        {
          for (int l=0; l<num_of_circles; l++)
          {
            double cy((dist+radius)+k*radius+2*k*dist);
            double cx((2+3*l)*radius);
            if  (sqrt((i-cx)*(i-cx)+(j-cy)*(j-cy)) < radius)
            {
              _M(i,j) = 1;
            }
          }
        }
      }
    }
    if (is_ok == false)
    {
      _M = _M.transpose();
    }

    std::cout << "Redistanciation en cours ... " << std::endl;
    LevelSet lv(_M);
    lv.redistancing(200);
    saveVTKFile(_M, nameMaskRedist);
  }
}
