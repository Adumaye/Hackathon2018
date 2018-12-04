// created by Gabriel Balestre, Stéphanie Lyon, Nathan Paillou, Théo Robushi
// supervised by Annabelle Collin and Heloise Beaugendre
// 2017/18

#ifndef UTIL_H
#define UTIL_H

#include "Dense"
#include <string>
#include <iostream>
#include <fstream>
#include <vector>

typedef Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic> field;

typedef struct
{
  double mu;
  double nu;
  double dt;
  double l1;
  double l2;
  double E;
  double neighbours;
  double coefFilter;
  std::string imageName;
  std::string scheme;
} config_t;


void parseFile(char* filename, config_t& c);

void saveVTKFile(const field& phi, std::string saveSolFileName);
void saveVTKFile(std::vector<std::vector<double>>& phi_v, std::string saveSolFileName);

std::vector<std::vector<double>> readVTKFile(std::vector<std::vector<double>>& phi_v, std::string fileName);

#endif
