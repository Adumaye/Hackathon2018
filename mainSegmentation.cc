// created by Gabriel Balestre, Stéphanie Lyon, Nathan Paillou, Théo Robuschi
// optimized by Baptiste Dufour, Antoine Dumaye, Marius Hérault, Odelin Gentieu et Louis Penin
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
#include "defclass.h"

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

  field u_0 = image->GetImage();

// Creation des myvector

  myvector<double> u0_myvector(rows*cols);

  for (size_t i = 0; i < rows; i++)
  {
    for (size_t j = 0; j < cols; j++)
    {
      u0_myvector[i*cols+j]= u_0(i,j);
    }
  }

  int nx=phi_v.size();
  int ny=phi_v[0].size();

  myvector<double> phi_myvector(nx*ny);

  for (int i=0; i<nx; i++)
  {
    for (int j=0; j<ny; j++)
    {
      phi_myvector[i*ny+j]= u_0(i,j);
    }
  }

  std::vector< std::vector <double> > newphi_v;

  newphi_v.resize(phi_v.size());
  for (int i=0;i<newphi_v.size() ;i++) { newphi_v[i].resize(phi_v[0].size()); }


  // newphi_myvector
  myvector<double> newphi_myvector(nx*ny);

  for (int i=0; i<nx; i++)
  {
    for (int j=0; j<ny; j++)
    {
      phi_myvector[i*ny+j]=phi_v[i][j];
    }
  }
// Fin de création de MyVector

  // Chan Vese method
  std::cout << "-------------------------------------------------" << std::endl;
  std::cout << "Initialisation de la méthode de Chan Vese" << std::endl;

  ChanVeseSchemes* chanVese=new ChanVeseSchemes(image);

  saveVTKFile(phi_v, "Results/sol_0.vtk");
  std::cout << "-------------------------------------------------" << std::endl;

  std::cout << "Iteration -- 1" << std::endl;

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
    newphi_myvector = chanVese->ExplicitScheme_myvector(phi_myvector,u0_myvector,c.dt,c.mu,c.nu,c.l1,c.l2, C1, C2, nx, ny);
    std::cout << "Fin myvector" << std::endl;

    std::cout << "taille de newphi_v" << newphi_v.size() << "et" << newphi_v[0].size() << std::endl;
    std::cout << "taille de newphi_myvector" << newphi_myvector.size() << std::endl;

    for (int i=0; i<nx; i++)
    {
      for (int j=0; j<ny;j++)
      {
        newphi_v[i][j]=newphi_myvector[i*ny+j];
      }
    }
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
    std::cout << "Fin CL" << std::endl;

    // Fin CL
  }
  else
  {
    std::cout << "Seulement le schéma ExplicitScheme est implémenté." << std::endl;
  }
  double diff = chanVese->fdiff(phi_v, newphi_v);

  phi_v = newphi_v;
  int i(2);
  // #pragma acc data copy(phi_v) create(newphi_v)
  {
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

        // #pragma acc data copy(newphi_v)
        // {
        //   #pragma acc parallel
        //   {
        //     #pragma acc loop
            for (int j=0; j<ny; ++j)
            {
              newphi_v[0][j]  = newphi_v[1][j];
              newphi_v[nx-1][j] = newphi_v[nx-2][j];
            }

            // #pragma acc loop
            for (int i=0; i<nx; ++i)
            {
              newphi_v[i][0]  = newphi_v[i][1];
              newphi_v[i][ny-1] = newphi_v[i][ny-2];
            }
        //   }
        //   #pragma acc update self(newphi_v)
        // }

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
        // #pragma acc update device(newphi_v)
        phi_v = newphi_v;
        i++;
      }
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


  std::cout<<" 0 "<< std::endl;
  u0_myvector.~myvector();
  std::cout<<" 1"<< std::endl;

  phi_myvector.~myvector();
  std::cout<<" 2 "<< std::endl;

  newphi_myvector.~myvector();
  std::cout<<" 3 "<< std::endl;


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
