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
  initMask->BuildMaskAndRedistancing(rows, cols, imagemaskdistance);
  std::cout << "-------------------------------------------------" << std::endl;

  // Chan Vese method
  std::cout << "-------------------------------------------------" << std::endl;
  std::cout << "Initialisation de la méthode de Chan Vese" << std::endl;
  ChanVeseSchemes* chanVese=new ChanVeseSchemes(image);
  field phi(initMask->GetRedistMask());
  saveVTKFile(phi, "Results/sol_0.vtk");
  std::cout << "-------------------------------------------------" << std::endl;

  std::cout << "Iteration -- 1" << std::endl;

  //Recopiage du phi de EIGEN dans un phi Vector
  std::vector< std::vector<double> >  phi_v;

  phi_v.resize(phi.rows());
  for (int i=0;i<phi.rows() ;i++) { phi_v[i].resize(phi.cols()); }

  for (int i=0 ; i < phi.rows(); i++)
  {
    for (int j=0 ; j < phi.cols(); j++)
    {
      phi_v[i][j]=phi(i,j);
    }
  }

  std::vector< std::vector <double> > newphi_v;

  newphi_v.resize(phi.rows());
  for (int i=0;i<newphi_v.size() ;i++) { phi_v[i].resize(phi.cols()); }

  std::string scheme(c.scheme);

  if (scheme == "ExplicitScheme")
  {
    std::cout << "Explicit scheme" << std::endl;
    newphi_v = chanVese->ExplicitScheme(phi_v,c.dt,c.mu,c.nu,c.l1,c.l2);

  }
  else
  {
    std::cout << "Seulement le schéma ExplicitScheme est implémenté." << std::endl;
  }
  double diff = chanVese->fdiff(phi_v, newphi_v);
  // double diff = (((newphi>=0).cast<double>()-0.5)*2. - ((phi>=0).cast<double>()-0.5)*2.).matrix().norm()
  // /(phi.rows()*phi.cols());
  phi_v = newphi_v;
  int i(2);
  while ( (diff > 5e-6) && (i < 100) )
  {
    if (i%10 == 0) { std::cout << "Iteration -- " << i << std::endl;}
    if (scheme == "ExplicitScheme")
    {
      newphi_v = chanVese->ExplicitScheme(phi_v,c.dt,c.mu,c.nu,c.l1,c.l2);

      diff = chanVese->fdiff(phi_v, newphi_v);
      // diff = (((newphi>=0).cast<double>()-0.5)*2. - ((phi>=0).cast<double>()-0.5)*2.).matrix().norm()
      // /(phi.rows()*phi.cols());
      if (i%10 == 0)
      {
        // newphi = ((newphi>=0).cast<double>()-0.5)*2;
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
  saveVTKFile(newphi, "Results/lastsol.vtk");
  image->WriteImage(chanVese->AbsGradPhi(newphi), (namewithoutextension + "_filtered_with_contour." + extension).c_str());

  cout << "Fin de la segmentation pour l'image." << endl;

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
