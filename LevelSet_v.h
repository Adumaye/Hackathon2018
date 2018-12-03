#ifndef LEVELSET_H_v
#define LEVELSET_H_v

//#include "Dense"
#include "Image.h"
#include "float.h"

// using namespace std;

class LevelSet_v
{
private:
    std::vector<std::vector<double>>& _phi_v;

    double sgn(double x) const;

    double pos_v(const std::vector<std::vector<double>>& X, int a, int b) const;
    double neg_v(const std::vector<std::vector<double>>& X, int a, int b) const;

    std::vector<std::vector<double>> godunov_v(const std::vector<std::vector<double>>& xm, const std::vector<std::vector<double>>& xp, const std::vector<std::vector<double>>& ym, const std::vector<std::vector<double>>& yp) const;

    double boundary(const double ij, const double imj,  const double ipj, const double ijm, const double ijp) const;

    void firstOrderUpwind_v(std::vector<std::vector<double>>& phi, const std::vector<std::vector<double>>& phi0, const double dt);

    double diff_m_v(const std::vector<std::vector<double>>& F, int index, int x, int y) const;

    double diff_p_v(const std::vector<std::vector<double>>& F, int index, int x, int y) const;

    std::vector<std::vector<double>> diff_m_v(const std::vector<std::vector<double>>& F, int index) const;

    std::vector<std::vector<double>> diff_p_v(const std::vector<std::vector<double>>& F, int index) const;

    std::vector<std::vector<double>> diff_v(const std::vector<std::vector<double>>&F, bool index) const;


  public:
      // LevelSet_v(std::vector<std::vector<double>>& phi_v);
      LevelSet_v(std::vector<std::vector<double>>& phi_v):_phi_v(phi_v){}
      ~LevelSet_v() = default;
      void redistancing_v(const int max_iter);
      void redistancing_v(const double eps);

      double fmax(const std::vector<std::vector<double>>& v) const;
      double fmin(const std::vector<std::vector<double>>& v) const;
};

#endif // LEVELSET_H_v
