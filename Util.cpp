// created by Gabriel Balestre, Stéphanie Lyon, Nathan Paillou, Théo Robushi
// supervised by Annabelle Collin and Heloise Beaugendre
// 2017/18

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>

#include "Util.h"

using namespace std;

// Lecture des fichiers de données
void parseFile(char* filename, config_t&c)
{
  ifstream file(filename);

  string line;
  while(getline(file, line))
  {
    stringstream iss(line);
    string key;
    if( getline(iss, key, '=') )
    {
      string value;
      getline(iss, value);
      const char * c_value = value.c_str();
      if(key == "mu")              c.mu = atof(c_value);
      else if(key == "nu")         c.nu = atof(c_value);
      else if(key == "dt")         c.dt = atof(c_value);
      else if(key == "l1")         c.l1 = atof(c_value);
      else if(key == "l2")         c.l2 = atof(c_value);
      else if(key == "name")       c.imageName = c_value;
      else if(key == "neighbours") c.neighbours = atoi(c_value);
      else if(key == "coefFilter") c.coefFilter = atof(c_value);
      else if(key == "scheme")     c.scheme = c_value;
    }
  }
}

// Écriture d'une matrice eigen ou d'un tableau dans un fichier
void saveVTKFile(const field& phi, std::string saveSolFileName)
{
  ofstream saveSol;
  saveSol.open(saveSolFileName, ios::out);

  saveSol << "# vtk DataFile Version 3.0" << endl;
  saveSol << "cell" << endl;
  saveSol << "ASCII" << endl;
  saveSol << "DATASET STRUCTURED_POINTS" << endl;
  saveSol << "DIMENSIONS " << phi.cols() << " " << phi.rows() << " " << 1 << endl;
  saveSol << "ORIGIN " << 0 << " " << 0 << " " << 0 << endl;
  saveSol << "SPACING " << 1.0 << " " << 1.0 << " " << 1 << endl;;
  saveSol << "POINT_DATA " <<  phi.cols()*phi.rows() << endl;
  saveSol << "SCALARS cell float" << endl;
  saveSol << "LOOKUP_TABLE default" << endl;

  for(int i=0; i<phi.rows(); i++)
  {
    for(int j=0; j<phi.cols(); j++)
    {
      saveSol << phi(i,j) << " ";
    }
  }
  saveSol << endl;
  saveSol.close();
}


// Écriture d'une matrice ou d'un tableau dans un fichier
void saveVTKFile(std::vector<std::vector<double>>& phi_v, std::string saveSolFileName)
{
  ofstream saveSol;
  saveSol.open(saveSolFileName, ios::out);

  saveSol << "# vtk DataFile Version 3.0" << endl;
  saveSol << "cell" << endl;
  saveSol << "ASCII" << endl;
  saveSol << "DATASET STRUCTURED_POINTS" << endl;
  saveSol << "DIMENSIONS " << phi_v.size() << " " << phi_v[0].size() << " " << 1 << endl;
  saveSol << "ORIGIN " << 0 << " " << 0 << " " << 0 << endl;
  saveSol << "SPACING " << 1.0 << " " << 1.0 << " " << 1 << endl;;
  saveSol << "POINT_DATA " <<  phi_v.size()*phi_v[0].size() << endl;
  saveSol << "SCALARS cell float" << endl;
  saveSol << "LOOKUP_TABLE default" << endl;

  for(int i=0; i<phi_v.size(); i++)
  {
    for(int j=0; j<phi_v[0].size(); j++)
    {
      saveSol << phi_v[i][j] << " ";
    }
  }
  saveSol << endl;
  saveSol.close();
}

void readVTKFile(field& phi, std::string fileName)
{
  ifstream file(fileName.data());
  if (!file.is_open())
  {
    cout << "Unable to open file " << fileName << endl;
    abort();
  }
  else
  {
    cout << "Reading file: " << fileName << endl;
  }
  string fileLine;
  while (!file.eof())
  {
    int Nx, Ny, Nz;
    std::string temp;
    getline(file, fileLine);

    if (fileLine.find("DATASET STRUCTURED_POINTS") != std::string::npos)
    {
      file >> temp >> Nx >> Ny >> Nz;
      phi.resize(Ny,Nx);
    }

    if (fileLine.find("LOOKUP_TABLE default") != std::string::npos)
    {
      for(int i=0; i<Ny; i++)
      {
        for(int j=0; j<Nx; j++)
        {
          file >> phi(i,j);
        }
      }
    }
  }
}
