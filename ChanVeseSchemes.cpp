// created by Gabriel Balestre, Stéphanie Lyon, Nathan Paillou, Théo Robushi
// supervised by Annabelle Collin and Heloise Beaugendre
// 2017/18

#ifndef _CHANVESESCHEMES_CPP
#define _CHANVESESCHEMES_CPP

#include "ChanVeseSchemes.h"
#include "float.h"
#include <iostream>
#include <fstream>
// #include <cmath>
// #include <mpi.h>
// #include <stdio.h>
// #include <stdlib.h>


ChanVeseSchemes::ChanVeseSchemes (Image* image) : _u0(image->GetImage())
{
	//création de _u0_v
	_u0_v.resize(_u0.rows());

	for (int i=0 ; i< _u0.rows(); i++)
	{
		_u0_v[i].resize(_u0.cols());
		for (int j=0 ; j< _u0.cols(); j++)
		{
			_u0_v[i][j]=_u0(i,j);
		}
	}
}

// Dérivées partielles
field ChanVeseSchemes::CSXPshift(const field& phi) const
{
	int nx(phi.rows()), ny(phi.cols());
	field CSXP(nx,ny); CSXP.leftCols(ny-1) = phi.rightCols(ny-1); CSXP.col(ny-1) = phi.col(0);
	return CSXP;
}
field ChanVeseSchemes::CSXMshift(const field& phi) const
{
	int nx(phi.rows()), ny(phi.cols());
	field CSXM(nx,ny); CSXM.rightCols(ny-1) = phi.leftCols(ny-1); CSXM.col(0) = phi.col(ny-1);
	return CSXM;
}
field ChanVeseSchemes::CSYPshift(const field& phi) const
{
	int nx(phi.rows()), ny(phi.cols());
	field CSYP(nx,ny); CSYP.topRows(nx-1) = phi.bottomRows(nx-1); CSYP.row(nx-1) = phi.row(0);
	return CSYP;
}
field ChanVeseSchemes::CSYMshift(const field& phi) const
{
	int nx(phi.rows()), ny(phi.cols());
	field CSYM(nx,ny); CSYM.bottomRows(nx-1) = phi.topRows(nx-1); CSYM.row(0) = phi.row(nx-1);
	return CSYM;
}
field ChanVeseSchemes::CSXMYPshift(const field& phi) const
{
	int nx(phi.rows()), ny(phi.cols());
	field CSXMYP(nx,ny); CSXMYP = CSYPshift(phi); CSXMYP = CSXMshift(CSXMYP);
	return CSXMYP;
}
field ChanVeseSchemes::CSXPYMshift(const field& phi) const
{
	int nx(phi.rows()), ny(phi.cols());
	field CSXPYM(nx,ny); CSXPYM = CSYMshift(phi); CSXPYM = CSXPshift(CSXPYM);
	return CSXPYM;
}
field ChanVeseSchemes::CSXMYMshift(const field& phi) const
{
	int nx(phi.rows()), ny(phi.cols());
	field CSXMYM(nx,ny); CSXMYM = CSYMshift(phi); CSXMYM = CSXMshift(CSXMYM);
	return CSXMYM;
}

// Discretisation du dirac
field ChanVeseSchemes::Dirac(const field& phi) const
{
	double eps(3.);
	return eps/(phi*phi+eps*eps);
}

// Discretisation de la valeur absolue du gradient
field ChanVeseSchemes::AbsGradPhi(const field& phi) const
{
	const double hx(1.), hy(1.0);

	field dxplus    = ( CSXPshift(phi) - phi ) / (hx);
	field dxminus   = ( phi - CSXMshift(phi) ) / (hx);
	field dyplus    = ( CSYPshift(phi) - phi ) / (hy);
	field dyminus   = ( phi - CSYMshift(phi) ) / (hy);
	field dxcentral = (dxplus+dxminus) / 2.;
	field dycentral = (dyplus+dyminus) / 2.;

	return sqrt(dxcentral*dxcentral + dycentral*dycentral);
}


double ChanVeseSchemes::fdiff(const std::vector<std::vector<double>>& phi_v, const std::vector<std::vector<double>>& newphi_v) const
{
	// diff = (((newphi>=0).cast<double>()-0.5)*2. - ((phi>=0).cast<double>()-0.5)*2.).matrix().norm()/(phi.rows()*phi.cols());
	double d=0.;

	for (int i=0; i<newphi_v.size(); i++)
	{
		for (int j=0; j<newphi_v[0].size(); j++)
		{
			d+= 4.*max( -phi_v[i][j]*newphi_v[i][j]/max(-phi_v[i][j]*newphi_v[i][j], 1.E-16) , 0. );
			//Petite modif par rapport au code d'origine : lorsque phi ou newphi était nul et l'autre étais négatif, d était incrémenté de 4.
		}
	}
	double diff = sqrt(d)/(newphi_v.size()*newphi_v[0].size());

	return diff;
}

std::pair<double,double> ChanVeseSchemes::Correction(const std::vector<std::vector<double>>& phi_v)
{
	// Calcul de C1 et C2
	std::pair<double,double> Correction;
	double dom_plus=0., dom_moins=0, z_plus=0., z_moins=0., C1, C2;

	for (int i=0; i<phi_v.size() ; i++)
	{
		for (int j=0; j<phi_v[0].size() ; j++)
		{
			dom_plus += max (phi_v[i][j]/max(abs(phi_v[i][j]),1.E-16),0.);
			dom_moins -= max (-phi_v[i][j]/max(abs(phi_v[i][j]),1.E-16),0.);
			//le max au dénominateur sert à ne jamais diviser par 0
			//dom_plus est le nombre d'éléments positifs dans phi
			//dom_moins est le nombre d'éléments négatifs dans phi
			z_plus += max (_u0_v[i][j]*phi_v[i][j]/max(abs(phi_v[i][j]),1.E-16),0.);
			z_moins -= max (-_u0_v[i][j]*phi_v[i][j]/max(abs(phi_v[i][j]),1.E-16),0.);
			//z_plus est la valeur de l'intégrale de z sur tous les phi(i,j) positifs
			//z_moins est la valeur de l'intégrale de z sur tous les phi(i,j) négatifs
		}
	}
	C1 = z_plus/dom_plus;
	C2 = z_moins/dom_moins;
	if (C1 < C2) // On sait pas trop pourquoi mais c'était fait comme ça dans la version précédente, donc on a fait pareil
	{
		double temp = C2;
		C2 = C1; C1 = temp;
	}
	// Fin de Calcul de C1 et C2
	Correction.first=C1;
	Correction.second=C2;
	return Correction;
}

std::vector<std::vector<double>>  ChanVeseSchemes::ExplicitScheme(const std::vector<std::vector<double>>&phi_v, const double dt,  const double mu, const double nu, const double l1, const double l2, const double C1, const double C2) const
{
  const double hx(1.), hy(1.0);
  const double eta(1e-8);

  int nx(phi_v.size());
  int ny(phi_v[0].size());

  std::vector< std::vector<double>> phiint_v;
  phiint_v.resize(nx);
  for (int i=0;i< nx;i++) { phiint_v[i].resize(ny); }

  double eps(3.);
  double diracij;
  //#pragma acc parallel
  {
    //#pragma acc loop independent
    for (int i=1; i<nx-1; ++i)
      {
	//#pragma acc loop independant
	for (int j=1; j<ny-1; ++j)
	  {
	    double firstterm   = (fdxplus(i,j,phi_v,hx)*coeffA(i,j,phi_v,hx,hy,eta) - fdxminus(i,j,phi_v,hx)*coeffA(i-1,j,phi_v,hx,hy,eta));
	    double secondterm  = (fdyplus(i,j,phi_v,hy)*coeffB(i,j,phi_v,hx,hy,eta) - fdyminus(i,j,phi_v,hy)*coeffB(i,j-1,phi_v,hx,hy,eta));
	    double correc      = -l1*(_u0(i,j)-C1)*(_u0(i,j)-C1) + l2*(_u0(i,j)-C2)*(_u0(i,j)-C2);
	    diracij            = eps/(pow(phi_v[i][j],2)+pow(eps,2));
	    phiint_v[i][j] = phi_v[i][j] + dt*diracij*(mu*(firstterm+secondterm)- nu + correc);
	  }
      }
  }

  return phiint_v;
}


#endif
