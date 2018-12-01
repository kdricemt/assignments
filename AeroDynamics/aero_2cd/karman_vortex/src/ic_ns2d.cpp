// main program
/******************************************************************************
                Incompressible Navier-Stokes 2D Flow Solver
                  "Flow Around Rectangular Cylinder"
                             K. Iiyama
                     (03-170313 U-Tokyo Aerospace)
******************************************************************************/

#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

#include "flowparam.hpp"
#include "setgrid.hpp"
#include "solveflow.hpp"
#include "streamline.hpp"

typedef std::vector<std::vector<double>> vector2double;
typedef std::vector<vector2double>       vector3double;

/***** Main *****/
int main() {
    std::ios_base::sync_with_stdio(false);

    /*int mdx = IcNs2d::SetGrid::mx + 100,
        mdy = IcNs2d::SetGrid::my + 100; */  // size of x[i][j], y[i][j]
    int           mx = SetGrid::mx;
    int           my = SetGrid::my;
    vector2double x(mx, std::vector<double>(my)), y(mx, std::vector<double>(my)),
        u(mx, std::vector<double>(my)), v(mx, std::vector<double>(my)),
        p(mx, std::vector<double>(my)), phi(mx, std::vector<double>(my));

    vector2double tx, ty;

    double dt = SetGrid::setgrd(x, y) * FlowParam::cfl;  // call setgrd(return min(dx, dy))

    int plot_step = 2000;  // 0: all steps until 5000

    std::cout << "*********************************************************** \n"
              << "                2D Incompressible Flow Solver \n"
              << "Reynolds Number " << FlowParam::re << "\n"
              << "Number of Grid Points: (x, y)= (" << SetGrid::mx << ", " << SetGrid::my << ")\n"
              << "CFL/dt/Steps= " << FlowParam::cfl << "/" << dt << "/" << FlowParam::nlast << "\n"
              << "Plot Step Number:" << plot_step << "\n"
              << "***********************************************************" << std::endl;

    SolveFlow::slvflw(dt, x, y, u, v, p, tx, ty, phi, plot_step);  // call slvflw

    int        nline_x0 = 90, nline_rear = 60, linelength = 250000;
    double     dtau = 0.0001;
    Streamline test(nline_x0, nline_rear, linelength, dtau);
    // int        dstep = 200, laststep = 5000;
    // test.streamline_all(dstep, laststep);  // for animation
    // int step = plot_step;
    // test.streamline_step(step);  // for figure image

    return 0;
}
/***** End Main *****/
