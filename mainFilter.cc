// created by Gabriel Balestre, Stéphanie Lyon, Nathan Paillou, Théo Robushi
// supervised by Annabelle Collin and Heloise Beaugendre
// 2017/18


#include <iostream>
#include <fstream>
#include <chrono>
#include <mpi.h>
#include "Image.h"
#include "InitMask.h"
#include "ChanVeseSchemes.h"
#include "Util.h"

using namespace std;

int main(int argc, char** argv)
{
  if (argc < 2)
  {
    cout << "Please, enter the name of your data file." << endl;
    abort();
  }
	config_t c;
	parseFile(argv[1],c);

	MPI_Init(NULL,NULL);
	int me, np;
	MPI_Comm_rank(MPI_COMM_WORLD, &me);
	MPI_Comm_size(MPI_COMM_WORLD, &np);

	Image* image=new Image(); //Pointeur vers l'image

  if (me == 0)
    cout << "Debut du calcul du filtre pour l'image : " << c.imageName << endl;

	image->ApplyMedianFilter(c.neighbours,c.coefFilter,c.imageName,me,np);

  if (me == 0)
    cout << "Fin du calcul du filtre pour l'image : " << c.imageName << endl;

  //image->WriteImage();

	return 0;
}
