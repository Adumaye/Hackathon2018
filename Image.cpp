// created by Gabriel Balestre, Stéphanie Lyon, Nathan Paillou, Théo Robushi
// supervised by Annabelle Collin and Heloise Beaugendre
// 2017/18

#ifndef _IMAGE_CPP

#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "mpi.h"
#include <tiffio.h>

#include "Image.h"

#undef Success
#include "Dense"

using namespace std;

Image::Image() {}

field Image::CutImage(int me, int np, std::string name)
{
  int Nx, Ny;
  Dimensions(name, Nx, Ny);
  int i0, i1, Nyloc;

  Load(me, i0, i1, Ny, np);
  Nyloc = i1-i0+1;
  _Im.resize(Nx,Nyloc);
  LocalLoading(_Im, i0, i1, name);
  return _Im;
}

void Image::Load(int me, int &i0, int &i1, int Ny, int np)
{
  int q = Ny/np;
  int r = Ny%np;
  if(r == 0){
    i0 = me*q;
    i1 = (me+1)*q-1;
  }
  else{
    if(me<r){
      i0 = me*q+me;
      i1 = q*(me+1)+me;
    }
    else{
      i0 = me*q+r;
      i1 = i0 + q-1;
    }
  }
}

int Image::GetI0(int me, int np, int Ny)
{
  int i0, i1;
  Load(me, i0, i1, Ny, np);
  return i0;
}

int Image::GetI1(int me, int np, int Ny)
{
  int i0, i1;
  Load(me, i0, i1, Ny, np);
  return i1;
}

void Image::Dimensions(std::string filename, int &Nx, int &Ny)
{
  TIFF *tif=TIFFOpen(filename.c_str(), "r");
  TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &Nx);
  TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &Ny);
}

void Image::BuildFilter(field phi, std::string filename, std::string output_file, int me, int np)
{
  int Nx, Ny;
  Dimensions(filename, Nx, Ny);
  TIFF * tif;

  bool test;
  test =false;
  int entree;
  entree = 0;

  while (test == false)
  {

    if (me == entree)
    {
      tif = TIFFOpen(output_file.c_str(), "a");
      TIFFSetDirectory(tif, 0);
      TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, Nx);
      TIFFSetField(tif, TIFFTAG_IMAGELENGTH, Ny);
      TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 8);
      TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
      TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
      TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 1);
      TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, 1);
      TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
      TIFFSetField(tif, TIFFTAG_RESOLUTIONUNIT, RESUNIT_NONE);
      TIFFSetField(tif, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
      if (me == 0)
      {
        unsigned char* v = new unsigned char[Nx];
        for(int j=0; j<Ny; j++)
        {
          for(int i=0; i<Nx; i++)
          {
            v[i] = 255;
          }
          TIFFWriteScanline(tif, v, j);
        }
      }
      int Nyloc, i0, i1;
      Load(me, i0, i1, Ny, np);
      Nyloc = i1-i0+1;

      unsigned char* v = new unsigned char[Nx];
      for(int j=0; j<Nyloc; j++)
      {
        for(int i=0; i<Nx; i++)
        {
          v[i] = phi(i,j);
        }
        TIFFWriteScanline(tif, v, j+i0);
      }
      TIFFClose(tif);
    }
    entree = entree + 1;
    if (entree==np)
    {
      test = true;
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
}


void Image::CreateFull(field phi, std::string filename, std::string output_file, int me, int np)
{
  int Nx, Ny,i0,i1, Nyloc;
  Dimensions(filename, Nx, Ny);
  TIFF * tif;

  bool test;
  test =false;
  int entree;
  entree = 0;

  Load(me,i0, i1,Ny,np);
  Nyloc =i1-i0+1;
  field M;
  M.resize(Nx, Nyloc);

  LocalLoading(M, i0, i1, filename);

  for (int j=1; j<Nyloc-1; j++)
  {
    for (int i=1; i<Nx-1; i++)
    {
      if( (phi(i,j)<0) && (phi(i-1,j)>0) )
      {M(i,j)=255;
      }
      else if ((phi(i,j)<0) && (phi(i+1,j)>0))
      {M(i,j)=255;
      }
      else if((phi(i,j)<0) && (phi(i,j-1)>0))
      {M(i,j)=255;
      }
      else if ((phi(i,j)<0) && (phi(i,j+1)>0))
      {M(i,j)=255;
      }
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);


  while (test == false)
  {

    if (me == entree)
    {
      tif = TIFFOpen(output_file.c_str(), "a");
      TIFFSetDirectory(tif, 0);
      TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, Nx);
      TIFFSetField(tif, TIFFTAG_IMAGELENGTH, Ny);
      TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 8);
      TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
      TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
      TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 1);
      TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, 1);
      TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
      TIFFSetField(tif, TIFFTAG_RESOLUTIONUNIT, RESUNIT_NONE);
      TIFFSetField(tif, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
      if (me == 0)
      {
        unsigned char* v = new unsigned char[Nx];
        for(int j=0; j<Ny; j++)
        {
          for(int i=0; i<Nx; i++)
          {
            v[i] = 255;
          }
          TIFFWriteScanline(tif, v, j);
        }
      }



      unsigned char* v = new unsigned char[Nx];
      for(int j=0; j<Nyloc; j++)
      {
        for(int i=0; i<Nx; i++)
        {
          v[i] = M(i,j);
        }
        TIFFWriteScanline(tif, v, j+i0);
      }
      TIFFClose(tif);
    }
    entree = entree + 1;
    if (entree==np)
    {
      test = true;
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
}

void Image::LocalLoading(field& M, int i0, int i1, string filename)
{
  int width, height, k, pourcentage;
  uint32* raster;
  const char *cstr = filename.c_str();

  TIFF *tif=TIFFOpen(cstr,"r");
  TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &width);           // uint32 width;
  TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &height);        // uint32 height;

  uint32 npixels=width*height;
  raster=(uint32 *) _TIFFmalloc(npixels *sizeof(uint32));
  TIFFReadRGBAImage(tif, width, height, raster, 0);
  k=0;
  pourcentage=0;

  for (int a =height-1 ; a!=-1 ; a=a-1)
  {
    for (int b=0 ; b!=width ; b=b+1)
    {
      if ((a>i0-1) && (a<i1+1))
      {
        uint16 R=(uint16 )TIFFGetB(raster[k]);
        M(b,a-i0)=R;
      }
      k+=1;
    }
  }
  _TIFFfree(raster);
  TIFFClose(tif);
}


void Image::ApplyMedianFilter(const int windows_size,const double seuil, string filename, int me, int np)
{
  MPI_Status Status2;
  bool test;
  test =false;

  int Nx, Ny, i0, i1, Nyloc, entree, tag;
  tag = 1000;
  entree = 0;
  Dimensions(filename, Nx, Ny);
  Load(me,i0, i1,Ny,np);
  Nyloc =i1-i0+1;

  double max, min;

  std::vector<double> windows_vector(int(pow(2*windows_size+1,2)));
  std::vector<double> windows_vector2;
  std::vector<double> moy_sig;

  // Matrice réduite de taille initiale
  field Im_reduit;
  Im_reduit.resize(Nx,Nyloc);

  _Im_u0.resize(Nx,Ny);




  //std::cout << "L'image traitée par le filtre médian est de taille " << Nx << "x" << Ny << std::endl;

  if ((me==0) || (me==np-1))
  {
    Nyloc = Nyloc+windows_size;
  }
  else
  {
    Nyloc = Nyloc+2*windows_size;
  }
  _Im_median.resize(Nx,Nyloc);
  _Im_origin.resize(Nx,Nyloc);

  if(me==0)
  {
    LocalLoading(_Im_origin, i0, i1+windows_size, filename);
  }
  else if (me==np-1)
  {
    LocalLoading(_Im_origin, i0-windows_size, i1, filename);
  }
  else
  {
    LocalLoading(_Im_origin, i0-windows_size, i1+windows_size, filename);
  }


  int pourcentage=0;

  int a,b;
  if (me==0)
  {
    a = 0;
    b = Nyloc-windows_size;
  }
  else if (me==np-1)
  {
    a = windows_size;
    b = Nyloc;
  }
  else
  {
    a = windows_size;
    b = Nyloc-windows_size;
  }

  for(int i=0; i<Nx; i++)
  {
    for(int j=a; j<b; j++)
    {
      if ((i<windows_size) || (i>=Nx-windows_size) || (j<windows_size) || (j>=Nyloc-windows_size))
      {
        _Im_median(i,j)=0;//_Im_origin(i,j); // c'est mieux de rien avoir que d'avoir les pixels d'origine...
      }
      else
      {
        for(int k=i-windows_size; k<i+windows_size+1; k++)
        {
          for(int l=j-windows_size; l<j+windows_size+1; l++)
          {
            windows_vector[k-i+windows_size+(l-j+windows_size)*(2*windows_size+1)]=_Im_origin(k,l);
          }
        }
        moy_sig=MeanAndStd(windows_vector);
        max=moy_sig[0]+seuil*moy_sig[1];
        min=moy_sig[0]-seuil*moy_sig[1];
        //si c'est du bruit
        if ((_Im_origin(i,j)>max) || (_Im_origin(i,j)<min))
        {
          windows_vector2.clear();
          for(int k=i-windows_size; k<i+windows_size+1; k++)
          {
            for(int l=j-windows_size; l<j+windows_size+1; l++)
            {
              //si autour c'est pas du bruit
              if ((_Im_origin(k,l)<=max) || (_Im_origin(k,l)>=min))
              {windows_vector2.push_back(_Im_origin(k,l));}
            }
          }
          _Im_median(i,j)=Median(windows_vector2);
        }
        else
        {
          _Im_median(i,j)=_Im_origin(i,j);
        }
      }
    }
    if (floor(i*100/Nx)>pourcentage)
    {
      pourcentage=pourcentage+1;
      if (me==0)
      {
	int i_barre;
	printf( "[" );
	for(i_barre=0;i_barre<=pourcentage;i_barre+=2) printf( "*" );
	for (;i_barre<=100; i_barre+=2 ) printf( "-" );
	printf( "] %3d %%", pourcentage );

	for(i_barre=0;i_barre<59;++i_barre) printf( "%c", 8 );

	fflush(stdout );
      }
    }
  }

  for(int i=0;i<Nx;i++)
  {
    for(int j=0; j<i1-i0+1; j++)
    {
      if (me==0)
      {
        Im_reduit(i,j)=_Im_median(i,j);
      }
      else
      {
        Im_reduit(i,j)=_Im_median(i,j+windows_size);
      }
      _Im_u0(i,j+i0) = Im_reduit(i,j);
    }
  }
  std::string extension(filename);
  extension = extension.substr(extension.find_last_of(".") + 1);
  std::string postprocessingimage(filename);
  postprocessingimage.erase(postprocessingimage.rfind('.'));
  postprocessingimage += "_filtered."+extension;
  BuildFilter(Im_reduit, filename, postprocessingimage, me, np);
}

void Image::ReadImage(std::string filename)
{
  TIFF* tif = TIFFOpen(filename.c_str(), "r");

  if (tif)
  {
    int width, height;
    uint32* raster;
    TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &width);           // uint32 width;
    TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &height);        // uint32 height;

    uint32 npixels=width*height;
    Eigen::VectorXd mat(npixels);

    raster=(uint32 *) _TIFFmalloc(npixels *sizeof(uint32));

    TIFFReadRGBAImage(tif, width, height, raster, 0);
    for (int k = 0 ; k < npixels ; k++)
    {
      uint16 R=(uint16) TIFFGetB(raster[k]);
      mat(k) = R;
    }
    _TIFFfree(raster);
    TIFFClose(tif);
    Eigen::Map<field> temp0((mat.segment(0,npixels)).data(),width,height);
    _Im_u0 = temp0;
  }
  else
  {
    std::cout << "Erreur lors de la lecture de l'image ! " << std::endl;
  }
}

void Image::WriteImage(const field& abs_grad_phi, std::string filename)
{
  TIFF* tif = TIFFOpen(filename.c_str(), "a");
  TIFFSetDirectory(tif, 0);
  TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, _Im_u0.rows());
  TIFFSetField(tif, TIFFTAG_IMAGELENGTH, _Im_u0.cols());
  TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 8);
  TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
  TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
  TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 1);
  TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, 1);
  TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
  TIFFSetField(tif, TIFFTAG_RESOLUTIONUNIT, RESUNIT_NONE);
  TIFFSetField(tif, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
  field full = 255*abs_grad_phi/abs_grad_phi.maxCoeff()+(1.-abs_grad_phi/abs_grad_phi.maxCoeff())*_Im_u0;

  for(int j=0; j<_Im_u0.cols(); j++)
  {
    unsigned char* v = new unsigned char[_Im_u0.rows()];
    for(int i=0; i<_Im_u0.rows(); i++)
    {
      v[i] = full(i,j);
    }
    TIFFWriteScanline(tif, v, j);
  }
  TIFFClose(tif);
}


std::vector<double> Image::MeanAndStd(const std::vector<double> & v) // Nathan méthode de Welford
{
  int i;
  int n(v.size());
  double x;
  double moyenne=0.;
  double old_moyenne=0;
  double sigma=0;
  for (i=0;i<n;i++)
  {
    x=v[i];
    old_moyenne=moyenne;
    moyenne=moyenne+(x-moyenne)/(i+1);
    sigma=sigma+(x-moyenne)*(x-old_moyenne);
  }
  sigma=sqrt(sigma/n);
  return {moyenne,sigma};
}

double Image::Median(std::vector<double> & v)
{
  //TRI PAR ORDRE CROISSANT
  double cur_min;//minimum courant
  int imin;//repère
  int i,j;//indice pour parcourir le tableau
  int n(v.size());
  for (i=0;i<n/2;i++)  // Nathan : rajouter /2
  {
    cur_min = v[i];
    imin = i;
    //repérer s'il existe un plus petit que lui dans le reste du tableau
    for (j=i+1;j<n;j++)
    {
      if (v[j]<cur_min)
      {
        imin = j;
        cur_min = v[j];
      }
    }
    v[imin] = v[i];
    v[i] = cur_min;
  }
  return v[n/2];
}

#define _IMAGE_CPP
#endif
