// created by Gabriel Balestre, Stéphanie Lyon, Nathan Paillou, Théo Robushi
// supervised by Annabelle Collin and Heloise Beaugendre
// 2017/18

#ifndef _IMAGE_H
#define _IMAGE_H

#include <string>
#include "Dense"
#include "Sparse"

typedef Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic> field;

class Image
{
private:
	field _Im; /*!< Black and white image of original image */
	field _Im_median; /*!< Filtered image with the median filter */
	field _Im_origin; /*!< Black and white image of original image */
	field _Im_u0;
public:

	Image();

  field CutImage (int me, int np, std::string name);

	void ReadImage(std::string filename);
	void WriteImage(const field& abs_grad_phi, std::string filename);
	void Load(int me, int &i0, int &i1, int Ny, int np);

	int GetI0(int me, int np, int Ny);
	int GetI1(int me, int np, int Ny);
	void Dimensions(std::string filename, int &Nx, int &Ny);
	void LocalLoading(field &M, int i0, int i1, std::string filename);
	void BuildFilter(field phi, std::string filename, std::string output_file, int me, int np);
	void CreateFull(field phi, std::string filename, std::string output_file, int me, int np);

	void ApplyMedianFilter(const int windows_size, const double seuil, std::string filename,int me,int np);
	double Median(std::vector<double> & v);
	std::vector<double> MeanAndStd(const std::vector<double> & v);

	const field & GetOriginalImage(){return _Im_origin;};
	const field & GetMedianImage(){return _Im_median;};
	const field & GetImage(){return _Im_u0;};
};

#endif
