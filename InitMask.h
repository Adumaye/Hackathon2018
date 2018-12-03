// created by Gabriel Balestre, Stéphanie Lyon, Nathan Paillou, Théo Robushi
// supervised by Annabelle Collin and Heloise Beaugendre
// 2017/18

#ifndef _MASK_H

#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

class InitMask
{
private:

public:
  InitMask();
  // Construit le masque et le redistancie si il n'existe pas déjà
  std::vector<std::vector<double>> BuildMaskAndRedistancing(int rows, int cols, std::string nameMaskRedist);
};

#define _MASK_H
#endif
