#ifndef SOLVEFLOW_H
#define SOLVEFLOW_H

#include "flowparam.hpp"
#include "setgrid.hpp"

typedef std::vector<std::vector<double>> vector2double;

// Solve Flow
class SolveFlow : private FlowParam, private SetGrid {
public:
    static void   slvflw(double &dt, vector2double &x, vector2double &y, vector2double &u,
                         vector2double &v, vector2double &p, vector2double &tx, vector2double &ty,
                         vector2double &phi, int plot_step);
    static void   intcnd(double &nbegin, double &time, vector2double &u, vector2double &v,
                         vector2double &p, vector2double &phi);
    static void   bcforp(vector2double &p);
    static void   bcforv(vector2double &u, vector2double &v);
    static void   bcforphi(vector2double &phi);
    static void   poiseq(double &dt, vector2double &u, vector2double &v, vector2double &p);
    static void   veloeq(double &dt, vector2double &u, vector2double &v, vector2double &p);
    static void   denseq(double &dt, vector2double &u, vector2double &v, vector2double &phi);
    static void   init_traj(vector2double &tx, vector2double &ty);
    static void   update_traj(vector2double &u, vector2double &v, vector2double &tx,
                              vector2double &ty, std::vector<int> &stop_cal_line, double &dt, int n);
    static double interp2(vector2double &v2, double x, double y, int &stop_cal);


private:
    /*static const int mx = SetGrid::mx, my = SetGrid::my, i1 = SetGrid::i1, i2 = SetGrid::i2,
                     j1 = SetGrid::j1, j2 = SetGrid::j2;  // copy values, dxi:1/dx*/
    static constexpr double dxi = 1.0 / dx, dyi = 1.0 / dy, dxi2 = dxi * dxi, dyi2 = dyi * dyi,
                            dxdy = 1.0 / (dxi * dxi + dyi * dyi) * 0.5,  //=2/(1/dx^2+1/dy^2)
        omegap = 1.1, errorp = 1e-4;
    static const int        maxitp   = 200;
    static const int        nline_x0 = 50, nline_side = 120;
    static constexpr double INF = 1000.0;
};

#endif
