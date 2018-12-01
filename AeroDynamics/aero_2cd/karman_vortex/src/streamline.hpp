/******************************************************************************
   program for calculating streamline
    filename for u,v datas: vmap10.txt,vmap20.txt,...
    files written in format like below
    x y u v
    x y u v
    x y u v
    .......
    x y u v

    x y u v <- new x
    x y u v
    .......
                                                        programmed by K.Iiyama
******************************************************************************/
#ifndef STREAMLINE_H
#define STREAMLINE_H

#include <string>
#include <vector>
#include "flowparam.hpp"
#include "setgrid.hpp"

typedef std::vector<std::vector<double>> vector2double;

class Streamline : private FlowParam, private SetGrid {
public:
    Streamline(int nline_x0, int nline_rear, int linelength, double dtau)
        : nline_(nline_x0 + nline_rear),
          nline_x0_(nline_x0),
          nline_rear_(nline_rear),
          linelength_(linelength),
          dtau_(dtau) {
        u_.resize(mx);
        v_.resize(mx);
        for (int i = 0; i < mx; i++) {
            u_[i].resize(my);
            v_[i].resize(my);
        }
        sx_.resize(nline_);
        sy_.resize(nline_);
        for (int i = 0; i < nline_; i++) {
            sx_[i].resize(linelength);
            sy_[i].resize(linelength);
        }
        std::cout << "Streamline plot" << std::endl;
        std::cout << "nline/nline_x0/nline_rear/linelength/dtau" << std::endl;
        std::cout << nline_ << " " << nline_x0 << " " << nline_rear << " " << linelength << " "
                  << dtau << std::endl;
        std::cout << "sxsize:" << sx_.size() << "x" << sx_[0].size() << std::endl;
        std::cout << "===================================================" << std::endl;
    };
    void streamline_all(int dstep, int laststep);
    void streamline_step(int step);

private:
    /*static const int        mx_ = SetGrid::mx, my_ = SetGrid::my;
    static const int        mdx_ = mx_ + 100, mdy_ = my_ + 100;
    static constexpr double dx_ = IcNs2d::SetGrid::dx, dy_ = IcNs2d::SetGrid::dy;
    static const int        i1_ = IcNs2d::SetGrid::i1, i2_ = IcNs2d::SetGrid::i2,
                     j1_ = IcNs2d::SetGrid::j1, j2_ = IcNs2d::SetGrid::j2;
    static const int  nlast_ = IcNs2d::nlast;*/
    const int         nline_;       // number of streamline
    const int         nline_x0_;    // number of streamline starting from x=0
    const int         nline_rear_;  // number of streamline starting from the rear of obstacle
    const int         linelength_;  // lenth of each streamline
    vector2double     u_, v_, sx_, sy_;
    const double      dtau_;
    const std::string file_dir_ = "./datafile/";

    void   readfile(int step);  // read file for steps
    void   update(int step);
    double interp2(vector2double& v2, double x, double y);
};
#endif
