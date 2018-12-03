// created by Gabriel Balestre, Stéphanie Lyon, Nathan Paillou, Théo Robushi
// supervised by Annabelle Collin and Heloise Beaugendre
// 2017/18

#ifndef _MASK_H

#include "Dense"
#include "LevelSet.h"

typedef Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic> field;

class InitMask
{
private:
  // Matrice représentant le masque redistancié
  field  _M ;

public:
  InitMask();
  // Construit le masque et le redistancie si il n'existe pas déjà
  void BuildMaskAndRedistancing(int rows, int cols, std::string nameMaskRedist);
  // Renvoie le masque
  const field & GetRedistMask(){return _M;};
};


#define _MASK_H
#endif
