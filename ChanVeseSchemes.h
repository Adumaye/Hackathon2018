// created by Gabriel Balestre, Stéphanie Lyon, Nathan Paillou, Théo Robuschi
// optimized by Baptiste Dufour, Antoine Dumaye, Marius Hérault, Odelin Gentieu et Louis Penin
// supervised by Annabelle Collin and Heloise Beaugendre
// 2017/18

#ifndef _CHANVESESCHEMES_H

#include "Dense"
#include "defclass.h"
#include "Image.h"
// #include "Util.h"
#include <vector>

using namespace std;

typedef Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic> field;

class ChanVeseSchemes
{


private:
	// Image à segmenter
	field _u0;
	std::vector<std::vector<double>> _u0_v;

public:
	ChanVeseSchemes(Image* image);

	// Dérivées partielles
	field CSXPshift(const field& phi) const;
	field CSXMshift(const field& phi) const;
	field CSYPshift(const field& phi) const;
	field CSYMshift(const field& phi) const;
	field CSXMYPshift(const field& phi) const;
	field CSXPYMshift(const field& phi) const;
	field CSXMYMshift(const field& phi) const;

	// Discretisation du Dirac
	field Dirac(const field& phi) const;

	// Valeur moyenne du domaine
	double ComputeMeanValueOnDomain(const field& phi) const;
	// Valeur moyenne sur le domaine complémentaire
	double ComputeMeanValueOnComplementaryDomain(const field& phi) const;

	// Valeur de diff -> condition d'arrêt de la boucle dans le mainSegmentation
	double fdiff(const std::vector<std::vector<double>> & phi_v, const std::vector<std::vector<double>>& newphi_v) const;

	// Correction
	std::pair<double,double> Correction(const std::vector<std::vector<double>>& phi_v);


	// |V phi|
	field AbsGradPhi(const field& phi) const;

	// Schéma pour différences finis
	// Explicit Scheme
	std::vector<std::vector<double>> ExplicitScheme(const std::vector<std::vector<double>> & phi_v, const double dt,  const double mu, const double nu,\
		 const double l1, const double l2, const double C1, const double C2) const;

	// myvector<double>  ExplicitScheme_myvector(myvector<double>& phi_v, myvector<double>& u0_myvector, const double dt,  const double mu, const double nu, const double l1,\
	// 	 const double l2, const double C1, const double C2, int nx, int ny) const;

	void ExplicitScheme_myvector(myvector<double>& newphi_v, const myvector<double>& phi_v, myvector<double>& u0_myvector, const double dt,  const double mu, const double nu, const double l1,\
		 const double l2, const double C1, const double C2, int nx, int ny) const;


	inline double fdxplus(int i,int j,const std::vector<std::vector<double>>& GrosPhi, double hx) const
	{
		return (GrosPhi[i+1][j]-GrosPhi[i][j])/hx;
	};

	inline double fdxminus(int i,int j, const std::vector<std::vector<double>>& GrosPhi, double hx) const
	{
		return (GrosPhi[i][j]-GrosPhi[i-1][j])/hx;
	};

	inline double fdyplus(int i,int j, const std::vector<std::vector<double>>& GrosPhi, double hy) const
	{
		return (GrosPhi[i][j+1]-GrosPhi[i][j])/hy;
	};

	inline double fdyminus(int i,int j, const std::vector<std::vector<double>>& GrosPhi, double hy) const
	{
		return (GrosPhi[i][j]-GrosPhi[i][j-1])/hy;
	};

	inline double fdxcentral(int i,int j, const std::vector<std::vector<double>>& GrosPhi, double hx) const
	{
		return (fdxplus(i,j,GrosPhi, hx)+fdxminus(i,j, GrosPhi,hx)) / 2.;
	};

	inline double fdycentral(int i,int j,const std::vector<std::vector<double>>& GrosPhi,double hy) const
	{
		return (fdyplus(i,j,GrosPhi, hy)+fdyminus(i,j,GrosPhi, hy)) / 2.;
	};

	inline double coeffA(int i,int j,const std::vector< std::vector<double>>& GrosPhi,double hx, double hy, const double eta) const
	{
		return 1./(sqrt(pow(eta,2) + pow(fdxplus(i,j,GrosPhi, hx),2) + pow(fdycentral(i,j,GrosPhi, hy),2)));
	};

	inline double coeffB(int i,int j,const std::vector< std::vector<double>>& GrosPhi,double hx, double hy, const double eta) const
	{
		return 1./(sqrt(pow(eta,2) + pow(fdyplus(i,j,GrosPhi, hy),2) + pow(fdxcentral(i,j,GrosPhi, hx),2)));
	};


	// //--------------DEBUT DES FONCTIONS MY_VECTOR-------------------------------------------------------------------------------

	inline double fdxplus_myvector(int pos,const myvector<double>& GrosPhi, double hx, int ny) const
	{
		return (GrosPhi[pos+ny]-GrosPhi[pos])/hx;
	};

	inline double fdxminus_myvector(int pos, const myvector<double>& GrosPhi, double hx, int ny) const
	{
		return (GrosPhi[pos]-GrosPhi[pos-ny])/hx;
	};

	inline double fdyplus_myvector(int pos,const myvector<double>& GrosPhi, double hy, int ny) const
	{
		return (GrosPhi[pos+1]-GrosPhi[pos])/hy;
	};

	inline double fdyminus_myvector(int pos,const myvector<double>& GrosPhi, double hy, int ny) const
	{
		return (GrosPhi[pos]-GrosPhi[pos-1])/hy;
	};

	inline double fdxcentral_myvector(int pos,const myvector<double>& GrosPhi, double hx, int ny) const
	{
		return (fdxplus_myvector(pos,GrosPhi, hx, ny)+fdxminus_myvector(pos,GrosPhi, hx, ny)) / 2.;
	};

	inline double fdycentral_myvector(int pos,const myvector<double>& GrosPhi,double hy, int ny) const
	{
		return (fdyplus_myvector(pos,GrosPhi, hy, ny)+fdyminus_myvector(pos, GrosPhi, hy, ny)) / 2.;
	};

	inline double coeffA_myvector(int pos,const myvector<double>& GrosPhi, double hx, double hy, const double eta, const double ny) const
	{
		return 1./(sqrt(pow(eta,2) + pow(fdxplus_myvector(pos, GrosPhi, hx, ny),2) + pow(fdycentral_myvector(pos, GrosPhi, hy, ny),2)));
	};

	inline double coeffB_myvector(int pos,const myvector<double>& GrosPhi, double hx, double hy, const double eta, const double ny) const
	{
		return 1./(sqrt(pow(eta,2) + pow(fdyplus_myvector(pos, GrosPhi, hy, ny),2) + pow(fdxcentral_myvector(pos, GrosPhi, hx, ny),2)));
	};


};

#define _CHANVESESCHEMES_H
#endif
