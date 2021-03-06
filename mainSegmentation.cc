// created by Gabriel Balestre, Stéphanie Lyon, Nathan Paillou, Théo Robushi
// supervised by Annabelle Collin and Heloise Beaugendre
// 2017/18


#include <iostream>
#include <fstream>
#include <chrono>
#include "Image.h"
#include "InitMask.h"
#include "ChanVeseSchemes.h"
#include "Util.h"
#include "LevelSet_v.h"

using namespace std;


int main(int argc, char** argv)
{
  // Lecture du fichier de données
  if (argc < 2)
  {
    cout << "Il faut indiquer le nom du fichier de données." << endl;
    abort();
  }
  config_t c;
  parseFile(argv[1],c);

  // Construction des différents noms de fichiers de sortie
  std::string extension(c.imageName);
  extension = extension.substr(extension.find_last_of(".") + 1);
  std::string namewithoutextension(c.imageName);
  namewithoutextension.erase(namewithoutextension.rfind('.'));
  c.imageName = namewithoutextension+"_filtered.tiff";
  std::string imagemaskdistance((namewithoutextension+"_filtered_distance_mask.vtk").c_str());

  // Créer un dossier
  system("mkdir -p ./Results");
  // et supprimer les anciens fichiers
  system("rm -f ./Results/*.vtk");

  // Lecture de l'image
  std::cout << "-------------------------------------------------" << std::endl;
  std::cout << "Lecture de l'image " << c.imageName << std::endl;
  Image* image = new Image();
  image->ReadImage(c.imageName);
  saveVTKFile(image->GetImage(), "Results/image.vtk");
  std::cout << "-------------------------------------------------" << std::endl;

  // Initialisation ou lecture du masque redistancié
  std::cout << "-------------------------------------------------" << std::endl;
  std::cout << "Construction ou lecture du masque redistancié" << std::endl;
  InitMask* initMask = new InitMask();
  int rows(image->GetImage().rows()), cols(image->GetImage().cols());
  std::vector< std::vector<double>> phi_v = initMask->BuildMaskAndRedistancing(rows, cols, imagemaskdistance);

  // Chan Vese method
  std::cout << "-------------------------------------------------" << std::endl;
  std::cout << "Initialisation de la méthode de Chan Vese" << std::endl;
  ChanVeseSchemes* chanVese=new ChanVeseSchemes(image);

  saveVTKFile(phi_v, "Results/sol_0.vtk");
  std::cout << "-------------------------------------------------" << std::endl;

  std::cout << "Iteration -- 1" << std::endl;

  std::vector< std::vector <double> > newphi_v;

  newphi_v.resize(phi_v.size());
  for (int i=0;i<newphi_v.size() ;i++) { phi_v[i].resize(phi_v[0].size()); }

  std::string scheme(c.scheme);

  double C1(0.);
  double C2(0.);
  std::pair<double,double> Correction;

//Début du Temps
  auto start = chrono::high_resolution_clock::now();

  if (scheme == "ExplicitScheme")
  {
    std::cout << "Explicit scheme" << std::endl;
    Correction = chanVese->Correction(phi_v);
    C1= Correction.first;
    C2= Correction.second;
    newphi_v = chanVese->ExplicitScheme(phi_v,c.dt,c.mu,c.nu,c.l1,c.l2, C1, C2);

    // CL
    int nx(phi_v.size());
  	int ny(phi_v[0].size());

    for (int j=0; j<ny; ++j)
    {
      newphi_v[0][j]  = newphi_v[1][j];
      newphi_v[nx-1][j] = newphi_v[nx-2][j];
    }

    for (int i=0; i<nx; ++i)
    {
      newphi_v[i][0]  = newphi_v[i][1];
      newphi_v[i][ny-1] = newphi_v[i][ny-2];
    }
    // Fin CL
  }
  else
  {
    std::cout << "Seulement le schéma ExplicitScheme est implémenté." << std::endl;
  }
  double diff = chanVese->fdiff(phi_v, newphi_v);

  phi_v = newphi_v;
  int i(2);
  while ( (diff > 5e-6) && (i < 100) )
  {
    if (i%10 == 0) { std::cout << "Iteration -- " << i << std::endl;}
    if (scheme == "ExplicitScheme")
    {
      Correction = chanVese->Correction(phi_v);
      C1= Correction.first;
      C2= Correction.second;
      newphi_v = chanVese->ExplicitScheme(phi_v,c.dt,c.mu,c.nu,c.l1,c.l2, C1, C2);

      // CL
      int nx(phi_v.size());
    	int ny(phi_v[0].size());

      for (int j=0; j<ny; ++j)
      {
        newphi_v[0][j]  = newphi_v[1][j];
        newphi_v[nx-1][j] = newphi_v[nx-2][j];
      }

      for (int i=0; i<nx; ++i)
      {
        newphi_v[i][0]  = newphi_v[i][1];
        newphi_v[i][ny-1] = newphi_v[i][ny-2];
      }
      // Fin CL

      diff = chanVese->fdiff(phi_v, newphi_v);

      if (i%10 == 0)
      {
        for (int i=0 ; i < newphi_v.size(); i++)
        {
          for (int j=0 ; j < newphi_v[0].size(); j++)
          {
            newphi_v[i][j]=newphi_v[i][j]/abs(newphi_v[i][j]);
          }
        }
        LevelSet_v lv_v(newphi_v);
        cout << "Redistanciation... " << endl;
        lv_v.redistancing_v(10);
        saveVTKFile(newphi_v, ("Results/sol_" + to_string(i)+ ".vtk").c_str());
        std::cout << "Evolution of phi : " << diff << std::endl;
      }
      phi_v = newphi_v;
      i++;
    }
  }

  field newphi(phi_v.size(),phi_v[0].size());
  for (int i=0 ; i < newphi_v.size(); i++)
  {
    for (int j=0 ; j < newphi_v[0].size(); j++)
    {
      newphi(i,j) = phi_v[i][j];
    }
  }

  newphi = ((newphi>=0).cast<double>()-0.5)*2;
  saveVTKFile(phi_v, "Results/lastsol.vtk");
  image->WriteImage(chanVese->AbsGradPhi(newphi), (namewithoutextension + "_filtered_with_contour." + extension).c_str());

  cout << "Fin de la segmentation pour l'image. !****Et bravo à Annabelle pour la petite Lise.****!" << endl;
  // Fin chrono
  auto finish = chrono::high_resolution_clock::now();
  double t = chrono::duration_cast<chrono::microseconds>(finish-start).count();
  cout << "Le prog a mis " << t*0.000001 << " secondes a s'effectuer" << endl;


  return 0;
}

// // CL
// int nx= phi_v.size();
// int ny= phi_v[0].size();
//
// for (int j=0; j<ny; ++j)
// {
//   newphi[0][j]  = newphi[1][j];
//   newphi[nx][j] = newphi[nx-1][j];
// }
//
// for (int i=0; i<nx; ++i)
// {
//   newphi[i][0]  = newphi[i][1];
//   newphi[i][ny] = newphi[i][ny-1];
// }
// // Fin CL
