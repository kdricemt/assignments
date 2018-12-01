#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include "streamline.hpp"

/***** readfile : read velocity data file for designated step *****/
void Streamline::readfile(int step) {
    std::string   filename = file_dir_ + "vmap/vmap" + std::to_string(step) + ".txt";
    std::ifstream ifs(filename);
    std::string   str;
    int           i = 0, j = 0;

    if (ifs.fail()) {
        std::cerr << "File do not exist.\n";
        exit(0);
    }

    while (getline(ifs, str)) {
        if (!str.empty()) {
            double x = 0, y = 0;
            std::sscanf(str.data(), "%lf %lf %lf %lf", &x, &y, &u_[i][j], &v_[i][j]);
            /*std::cout << step << " " << x << " " << y << " " << u_[i][j] << " " << v_[i][j] << " "
                      << std::endl;*/
            j += 1;
        } else {
            i += 1;
            j = 0;
        }
    }
}
/***** end readfile *****/

void Streamline::update(int step) {
    std::cout << "calclating step:" << step << std::endl;
    readfile(step);
    // initialize lines startting from x=0
    std::cout << "x=0" << std::endl;
    for (int line = 0; line < nline_x0_; line++) {
        double y_scale = 0.2;
        sx_[line][0]   = -(i1 + i2) * dx / 2;
        sy_[line][0]   = (dy * my * ((double)line / nline_x0_) - (j1 + j2) * dy / 2) * y_scale;
    }
    // initialize lines starting from the rear of the obstacle, starts from 3 point groups
    int nline_rear1 = nline_rear_ / 3;  // number of points in the line
    int nline_rear2 = nline_rear_ / 3;
    int nline_rear3 = nline_rear_ - nline_rear1 - nline_rear2;

    /*     |->  |->  |->
           |->  |->  |->
           |->  |->  |->
           |->  |->  |->
           |->  |=>  |->   3 start point lines behind the obstacle   */

    double x_offset_1 = 0.01;  // distance from the rear edge of the obstacle
    double y_scale_1  = 1.5;   // length of the line / ((j2-j1)*dy)
    double y_offset_1 = 0;     // center of the line (y)

    double x_offset_2 = 0.5;
    double y_scale_2  = 1.5;
    double y_offset_2 = 0;

    double x_offset_3 = 1;
    double y_scale_3  = 1.5;
    double y_offset_3 = 0;

    std::cout << "Rear1" << std::endl;
    for (int line = nline_x0_; line < nline_x0_ + nline_rear1; line++) {
        sx_[line][0] = ((i2 - i1) / 2) * dx + x_offset_1;
        sy_[line][0] = y_offset_1 + ((j2 - j1) * ((double)line - (double)nline_x0_) / nline_rear1 -
                                     ((j1 + j2) / 2 - j1)) *
                                        dy * y_scale_1;
    }

    std::cout << "Rear2" << std::endl;
    for (int line = nline_x0_ + nline_rear1; line < nline_x0_ + nline_rear1 + nline_rear2; line++) {
        sx_[line][0] = ((i2 - i1) / 2) * dx + x_offset_2;
        sy_[line][0] =
            y_offset_2 +
            ((j2 - j1) * ((double)line - (double)nline_x0_ - (double)nline_rear1) / nline_rear2 -
             ((j1 + j2) / 2 - j1)) *
                dy * y_scale_2;
    }

    std::cout << "Rear3" << std::endl;
    for (int line = nline_x0_ + nline_rear1 + nline_rear2; line < nline_; line++) {
        sx_[line][0] = ((i2 - i1) / 2) * dx + x_offset_3;
        sy_[line][0] =
            y_offset_3 +
            ((j2 - j1) *
                 ((double)line - (double)nline_x0_ - (double)nline_rear1 - (double)nline_rear2) /
                 nline_rear3 -
             ((j1 + j2) / 2 - j1)) *
                dy * y_scale_3;
    }

    // runge-kutta
    for (int line = 0; line < nline_; line++) {
        for (int m = 0; m < linelength_ - 1; m++) {
            double xm = sx_[line][m], ym = sy_[line][m];
            /*std::cout << "(step,line,m)=" << step << "," << line << "," << m << "    (xm,ym)=" <<
               xm
                      << " " << ym << std::endl;*/

            double kx1 = dtau_ * interp2(u_, xm, ym);
            double ky1 = dtau_ * interp2(v_, xm, ym);

            double kx2 = dtau_ * interp2(u_, xm + kx1 / 2, ym + ky1 / 2);
            double ky2 = dtau_ * interp2(v_, xm + kx1 / 2, ym + ky1 / 2);

            double kx3 = dtau_ * interp2(u_, xm + kx2 / 2, ym + ky2 / 2);
            double ky3 = dtau_ * interp2(v_, xm + kx2 / 2, ym + ky2 / 2);

            double kx4 = dtau_ * interp2(u_, xm + kx3, ym + ky3);
            double ky4 = dtau_ * interp2(v_, xm + kx3, ym + ky3);

            sx_[line][m + 1] = xm + (kx1 + 2 * kx2 + 2 * kx3 + kx4) / 6;
            sy_[line][m + 1] = ym + (ky1 + 2 * ky2 + 2 * ky3 + ky4) / 6;
        }
    }
}

double Streamline::interp2(vector2double& v2, double x, double y) {
    int i0 = (i1 + i2) / 2, j0 = (j1 + j2) / 2;  // index of (0,0)
    int xi = (int)(x / dx), yi = (int)(y / dy);  // index offset
    if (x < 0 && x > -10) xi += -1;
    if (y < 0 && y > -10) yi += -1;
    int    i = xi + i0, j = yi + j0;  // index of (x,y)
    double s      = (x - xi * dx) / dx,
           t      = (y - yi * dy) / dy;  // coordinate/coordnate = []
    double bottom = (1 - s) * v2[i][j] + s * v2[i + 1][j];
    double top    = (1 - s) * v2[i][j + 1] + s * v2[i + 1][j + 1];
    // std::cout << "x/y/xi/yi/s/t" << std::endl;
    // std::cout << x << " " << y << " " << i << " " << j << " " << s << " " << t << std::endl;
    return (1 - t) * bottom + t * top;
};

void Streamline::streamline_all(int dstep, int laststep) {
    if (laststep > nlast || (dstep % 10)) {
        std::cerr << "reenter dstep or laststep.\n";
        exit(0);
    } else {
        for (int step = 10; step < laststep; step = step + dstep) {
            update(step);
            for (int line = 0; line < nline_; line++) {
                std::fstream fout;
                std::string  filename = file_dir_ + "streamline/" + std::to_string(step) + "_" +
                                       std::to_string(line) + ".txt";
                fout.open(filename, std::ios::out);

                for (int length = 0; length < linelength_; length++) {
                    fout << sx_[line][length] << " " << sy_[line][length] << std::endl;
                }
                fout.close();
            }
        }
    }
}

void Streamline::streamline_step(int step) {
    if (step % nlp) {
        std::cerr << "reenter step.\n";
        exit(0);
    } else {
        update(step);
        for (int line = 0; line < nline_; line++) {
            std::fstream fout;
            std::string  filename = file_dir_ + "streamline/" + std::to_string(step) + "_" +
                                   std::to_string(line) + ".txt";
            fout.open(filename, std::ios::out);
            std::cout << "line number:" << line << std::endl;

            for (int length = 0; length < linelength_; length++) {
                fout << sx_[line][length] << " " << sy_[line][length] << std::endl;
            }
            fout.close();
        }
    }
}
