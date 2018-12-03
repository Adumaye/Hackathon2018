#ifndef LEVELSET_H
#define LEVELSET_H

//#include "Dense"
#include "Image.h"
#include <iostream>


typedef Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic> field;

class LevelSet
{
private:
    field& _phi;

    std::vector<field> normal() const;
    field curvature(std::vector<field> &n) const;

    double sgn(double x) const;
    field sgn_peng() const;
    field sgn_sussman(const field &X) const;

    field pos(const field& X) const;
    field neg(const field& X) const;

    field godunov(const field& xm, const field& xp, const field& ym, const field& yp) const;

    double boundary(const double ij, const double imj,  const double ipj, const double ijm, const double ijp) const;

    void firstOrderUpwind(field& phi, const field& phi0, const double dt);


    double diff_m(const field& F, int index, int x, int y) const;


    double diff_p(const field& F, int index, int x, int y) const;

    field diff_m(const field& F, int index) const;

    field diff_p(const field& F, int index) const;

    field diff(const field &F, bool index) const;


  public:
      LevelSet(field& phi):_phi(phi){}
      ~LevelSet() = default;
      void redistancing(const int max_iter);
      void redistancing(const double eps);

      field curvature() const;
      field check() const;
      bool check(const double eps) const;

};

#endif // LEVELSET_H
