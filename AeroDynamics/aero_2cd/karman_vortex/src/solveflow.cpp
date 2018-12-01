#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include "setgrid.hpp"
#include "solveflow.hpp"


/*** slvflw : write pv plot data and trajectory data
when plot_step = 0, write down pv plot data for all steps  ***/
void SolveFlow::slvflw(double &dt, vector2double &x, vector2double &y, vector2double &u,
                       vector2double &v, vector2double &p, vector2double &tx, vector2double &ty,
                       vector2double &phi, int plot_step) {
    if ((plot_step % nlp)) {
        std::cerr << "reenter step.\n";
        exit(0);
    }
    double           nbegin, time;
    std::vector<int> stop_cal_line(nline_x0 + nline_side);

    for (int i = 0; i < nline_x0 + nline_side; i++) {
        stop_cal_line[i] = 0;  // array for flags; 1 if the line reached the right edge
    }

    std::cout << "solveflow" << std::endl;
    // set initial conditions
    intcnd(nbegin, time, u, v, p, phi);
    bcforp(p);
    bcforv(u, v);  // set boundary conditions for initial conditions
    bcforphi(phi);
    init_traj(tx, ty);

    //  cout<<" "<<dxi<<" "<<dxi2<<endl;
    std::cout << "Step/cd/cl/cp1/cp2" << std::endl;

    for (int n = 0; n <= nlast; n++) {
        double nstep = n + nbegin;
        time += dt;

        denseq(dt, u, v, phi);                            // update density
        bcforphi(phi);                                    // set boundary conditions for density
        poiseq(dt, u, v, p);                              // update pressure
        bcforp(p);                                        // set boundary conditions for pressure
        veloeq(dt, u, v, p);                              // update velocity
        bcforv(u, v);                                     // set boundary conditions for velocity
        update_traj(u, v, tx, ty, stop_cal_line, dt, n);  // update trajectory

        // calculate CL and CD
        double cd = 0.0;
        for (int j = j1; j <= j2 - 1; j++) {
            double cpfore = 1.0 * (p[i1][j] + p[i1][j + 1]);  // integrate cp
            double cpback = 1.0 * (p[i2][j] + p[i2][j + 1]);
            cd += SetGrid::dy * (cpfore - cpback);
        }

        double cl = 0.0;
        for (int i = i1; i <= i2 - 1; i++) {
            double cpbtm = 1.0 * (p[i][j1] + p[i + 1][j1]);
            double cptop = 1.0 * (p[i][j2] + p[i + 1][j2]);
            cl += SetGrid::dx * (cpbtm - cptop);
        }

        // display and write results for each nlp steps
        if (!(n % nlp)) {
            double cp1 = 2.0 * p[i2 + i2 - i1][j1];
            double cp2 = 2.0 * p[i2 + i2 - i1][j2];
            std::cout << std::fixed << std::setprecision(2) << nstep << " " << cd << " " << cl
                      << " " << cp1 << " " << cp2 << std::endl;
            // data file for pressure
            std::fstream fout_p;
            std::string  filename_p = "./datafile/pmap/pmap" + std::to_string(n) + ".txt";
            fout_p.open(filename_p, std::ios::out);
            // data file for velocity
            std::fstream fout_v;
            std::string  filename_v = "./datafile/vmap/vmap" + std::to_string(n) + ".txt";
            fout_v.open(filename_v, std::ios::out);
            // data file for density
            std::fstream fout_d;
            std::string  filename_d = "./datafile/density/density" + std::to_string(n) + ".txt";
            fout_d.open(filename_d, std::ios::out);
            // write data file for Trajectory
            for (int line = 0; line < nline_x0 + nline_side; line++) {
                std::fstream fout_tr;
                std::string  filename_tr =
                    "./datafile/trajectory/trajectory" + std::to_string(line) + ".txt";
                if (n == 0)
                    fout_tr.open(filename_tr, std::ios::out);
                else
                    fout_tr.open(filename_tr, std::ios::app);
                fout_tr << tx[line][n] << " " << ty[line][n] << std::endl;
                fout_tr.close();
            }
            // write results for all steps
            if (plot_step == 0) {
                for (int i = 0; i < mx; i++) {
                    for (int j = 0; j < my; j++) {
                        fout_p << x[i][j] << " " << y[i][j] << " " << 2.0 * p[i][j] << std::endl;
                        fout_v << x[i][j] << " " << y[i][j] << " " << u[i][j] << " " << v[i][j]
                               << std::endl;
                        fout_d << x[i][j] << " " << y[i][j] << " " << phi[i][j] << std::endl;
                    }
                    fout_p << std::endl;
                    fout_v << std::endl;
                    fout_d << std::endl;
                }
                fout_p.close();
                fout_v.close();
                fout_d.close();
            }
            // write only until step = plot_step
            else if (n <= plot_step) {
                for (int i = 0; i < mx; i++) {
                    for (int j = 0; j < my; j++) {
                        fout_p << x[i][j] << " " << y[i][j] << " " << 2.0 * p[i][j] << std::endl;
                        fout_v << x[i][j] << " " << y[i][j] << " " << u[i][j] << " " << v[i][j]
                               << std::endl;
                        fout_d << x[i][j] << " " << y[i][j] << " " << phi[i][j] << std::endl;
                    }
                    fout_p << std::endl;
                    fout_v << std::endl;
                    fout_d << std::endl;
                }
                fout_p.close();
                fout_v.close();
                fout_d.close();

                if (n == plot_step) {
                    std::cout << "finish calculation -- Step No: " << n << std::endl;
                    break;
                }
            }
        }
    }
}
/***** End slvflw *****/

/*** inittraj : trajectory lines starting point ***/
void SolveFlow::init_traj(vector2double &tx, vector2double &ty) {
    tx.resize(nline_x0 + nline_side);
    ty.resize(nline_x0 + nline_side);
    for (int i = 0; i < nline_x0 + nline_side; i++) {
        tx[i].resize(nlast + 2);
        ty[i].resize(nlast + 2);
    }
    std::cout << tx.size() << " " << tx[0].size() << std::endl;
    std::cout << "init x=0" << std::endl;
    for (int line = 0; line < nline_x0; line++) {
        double y_scale = 0.2;
        tx[line][0]    = -(i1 + i2) * dx / 2;
        ty[line][0]    = (dy * my * ((double)line / nline_x0) - (j1 + j2) * dy / 2) * y_scale;
    }
    // initialize lines starting from the surrounding points of the obstacle
    double offset      = 0.01;  // distance from the sides of the obstacle
    double rec_sidelen = (i2 - i1) * dx;
    int    index_side1 = nline_x0;
    int    index_side2 = nline_x0 + nline_side / 4;
    int    index_side3 = nline_x0 + nline_side * 2 / 4;
    int    index_side4 = nline_x0 + nline_side * 3 / 4;

    std::cout << "init  side" << std::endl;
    for (int line = index_side1; line < index_side2; line++) {
        tx[line][0] = -rec_sidelen / 2 - offset;
        ty[line][0] =
            rec_sidelen / 2 + offset - (rec_sidelen / (nline_side / 4)) * (line - index_side1);
    }
    for (int line = index_side2; line < index_side3; line++) {
        tx[line][0] =
            -rec_sidelen / 2 - offset + (rec_sidelen / (nline_side / 4)) * (line - index_side2);
        ty[line][0] = -rec_sidelen / 2 - offset;
    }
    for (int line = index_side3; line < index_side4; line++) {
        tx[line][0] = rec_sidelen / 2 + offset;
        ty[line][0] =
            -rec_sidelen / 2 - offset + (rec_sidelen / (nline_side / 4)) * (line - index_side3);
    }
    for (int line = index_side4; line < nline_x0 + nline_side; line++) {
        tx[line][0] =
            rec_sidelen / 2 + offset - (rec_sidelen / (nline_side / 4)) * (line - index_side4);
        ty[line][0] = rec_sidelen / 2 + offset;
    }
}


void SolveFlow::update_traj(vector2double &u, vector2double &v, vector2double &tx,
                            vector2double &ty, std::vector<int> &stop_cal_line, double &dt, int n) {
    // runge-kutta
    for (int line = 0; line < nline_x0 + nline_side; line++) {
        int    stop_cal     = 0;
        int    already_stop = 0;
        double kx1 = 0, ky1 = 0, kx2 = 0, ky2 = 0, kx3 = 0, ky3 = 0, kx4 = 0, ky4 = 0;
        double xn = 0, yn = 0;
        if (stop_cal_line[line] == 1) {
            already_stop = 1;
        }
        if (already_stop == 0) {
            xn = tx[line][n];
            yn = ty[line][n];
            /*std::cout << "(step,line,m)=" << step << "," << line << "," << m << "    (xm,ym)=" <<
               xm        << " " << ym << std::endl;*/

            kx1 = dt * interp2(u, xn, yn, stop_cal);
            ky1 = dt * interp2(v, xn, yn, stop_cal);

            kx2 = dt * interp2(u, xn + kx1 / 2, yn + ky1 / 2, stop_cal);
            ky2 = dt * interp2(v, xn + kx1 / 2, yn + ky1 / 2, stop_cal);

            kx3 = dt * interp2(u, xn + kx2 / 2, yn + ky2 / 2, stop_cal);
            ky3 = dt * interp2(v, xn + kx2 / 2, yn + ky2 / 2, stop_cal);

            kx4 = dt * interp2(u, xn + kx3, yn + ky3, stop_cal);
            ky4 = dt * interp2(v, xn + kx3, yn + ky3, stop_cal);
        }
        if (stop_cal == 1 || already_stop == 1) {
            if (already_stop == 0) {
                std::cout << "Line No: " << line << " --reached the right edge" << std::endl;
                stop_cal_line[line] = 1;
            }
            tx[line][n + 1] = tx[line][n];
            ty[line][n + 1] = ty[line][n];
        } else {
            tx[line][n + 1] = xn + (kx1 + 2 * kx2 + 2 * kx3 + kx4) / 6;
            ty[line][n + 1] = yn + (ky1 + 2 * ky2 + 2 * ky3 + ky4) / 6;
        }
    }
}

double SolveFlow::interp2(vector2double &v2, double x, double y, int &stop_cal) {
    int i0 = (i1 + i2) / 2, j0 = (j1 + j2) / 2;  // index of (0,0)
    int xi = (int)(x / dx), yi = (int)(y / dy);  // index offset
    if (x < 0 && x > -10) xi += -1;
    if (y < 0 && y > -10) yi += -1;
    int i = xi + i0, j = yi + j0;  // index of (x,y)
    if (i >= mx - 1) {
        stop_cal = 1;
        return INF;
    }
    double s      = (x - xi * dx) / dx,
           t      = (y - yi * dy) / dy;  // coordinate/coordnate = []
    double bottom = (1 - s) * v2[i][j] + s * v2[i + 1][j];
    double top    = (1 - s) * v2[i][j + 1] + s * v2[i + 1][j + 1];
    // std::cout << "x/y/xi/yi/s/t" << std::endl;
    // std::cout << x << " " << y << " " << i << " " << j << " " << s << " " << t << std::endl;
    return (1 - t) * bottom + t * top;
};


/***** intcnd : set initial conditions *****/
void SolveFlow::intcnd(double &nbegin, double &time, vector2double &u, vector2double &v,
                       vector2double &p, vector2double &phi) {
    nbegin = 0;
    time   = 0.0;

    for (int i = 0; i < mx; i++) {
        for (int j = 0; j < my; j++) {
            u[i][j]   = 1.0;
            v[i][j]   = 0.0;
            p[i][j]   = 0.0;
            phi[i][j] = 0.0;
        }
    }

    std::cout << "intcnd" << std::endl;
}
/***** End intcnd *****/

/***** bcforp: set boundary conditions for pressure *****/
void SolveFlow::bcforp(vector2double &p) {
    //  cout<<" bcforphoge";
    for (int j = 0; j < my; j++) {
        int i   = 1;
        p[i][j] = 0.0;  // inflow(i = 0)
        i       = mx - 1;
        p[i][j] = 0.0;  // downstream (i = mx - 1)
    }

    for (int i = 0; i < mx; i++) {
        int j   = 1;
        p[i][j] = 0.0;  // bottom (j = 0)
        j       = my - 1;
        p[i][j] = 0.0;  // top (j = my - 1)
    }

    p[i1][j1] = p[i1 - 1][j1 - 1];  // wall condition
    p[i1][j2] = p[i1 - 1][j2 + 1];
    p[i2][j1] = p[i2 + 1][j1 - 1];
    p[i2][j2] = p[i2 + 1][j2 + 1];

    for (int j = j1 + 1; j <= j2 - 1; j++) {
        p[i1][j] = p[i1 - 1][j];
        p[i2][j] = p[i2 + 1][j];
    }
    for (int i = i1 + 1; i <= i2 - 1; i++) {
        p[i][j1] = p[i][j1 - 1];
        p[i][j2] = p[i][j2 + 1];
    }
}
/***** End bcforp *****/

/***** bcforv: set boundary conditions for velocity *****/
void SolveFlow::bcforv(vector2double &u, vector2double &v) {
    for (int j = 0; j < my; j++) {
        int i       = 1;
        u[i][j]     = 1;  // inflow
        v[i][j]     = 0.0;
        u[i - 1][j] = 1;
        v[i - 1][j] = 0.0;

        i           = mx - 2;
        u[i][j]     = 2.0 * u[i - 1][j] - u[i - 2][j];  // downstream
        v[i][j]     = 2.0 * v[i - 1][j] - v[i - 2][j];
        u[i + 1][j] = 2.0 * u[i][j] - u[i - 1][j];
        v[i + 1][j] = 2.0 * v[i][j] - v[i - 1][j];
    }

    for (int i = 0; i < mx; i++) {
        int j       = 1;
        u[i][j]     = 2.0 * u[i][j + 1] - u[i][j + 2];  // bottom
        v[i][j]     = 2.0 * v[i][j + 1] - v[i][j + 2];
        u[i][j - 1] = 2.0 * u[i][j] - u[i][j + 1];
        v[i][j - 1] = 2.0 * v[i][j] - v[i][j + 1];

        j           = my - 2;
        u[i][j]     = 2.0 * u[i][j - 1] - u[i][j - 2];  // top
        v[i][j]     = 2.0 * v[i][j - 1] - v[i][j - 2];
        u[i][j + 1] = 2.0 * u[i][j] - u[i][j - 1];
        v[i][j + 1] = 2.0 * v[i][j] - v[i][j - 1];
    }

    for (int i = i1; i <= i2; i++) {  // wall condition
        for (int j = j1; j <= j2; j++) {
            u[i][j] = 0.0;
            v[i][j] = 0.0;
        }
    }
}
/***** End bcforv *****/

/***** bcforp: set boundary conditions for density *****/
void SolveFlow::bcforphi(vector2double &phi) {
    //  cout<<" bcforphoge";
    for (int j = 0; j < my; j++) {
        int i = 1;
        if ((j > j1 + 1 && j < j1 + 4) || (j > j2 - 4 && j < j2 - 1)) {
            phi[i][j]     = 0.5;  // inflow
            phi[i - 1][j] = 0.5;
        }

        i             = mx - 2;
        phi[i][j]     = 2.0 * phi[i - 1][j] - phi[i - 2][j];  // downstream
        phi[i + 1][j] = 2.0 * phi[i][j] - phi[i - 1][j];
    }

    for (int i = 0; i < mx; i++) {
        int j         = 1;
        phi[i][j]     = 2.0 * phi[i][j + 1] - phi[i][j + 2];  // bottom
        phi[i][j - 1] = 2.0 * phi[i][j] - phi[i][j + 1];

        j             = my - 2;
        phi[i][j]     = 2.0 * phi[i][j - 1] - phi[i][j - 2];  // top
        phi[i][j + 1] = 2.0 * phi[i][j] - phi[i][j - 1];
    }

    phi[i1][j1] = phi[i1 - 1][j1 - 1];  // wall condition
    phi[i1][j2] = phi[i1 - 1][j2 + 1];
    phi[i2][j1] = phi[i2 + 1][j1 - 1];
    phi[i2][j2] = phi[i2 + 1][j2 + 1];

    for (int j = j1 + 1; j <= j2 - 1; j++) {
        phi[i1][j] = phi[i1 - 1][j];
        phi[i2][j] = phi[i2 + 1][j];
    }
    for (int i = i1 + 1; i <= i2 - 1; i++) {
        phi[i][j1] = phi[i][j1 - 1];
        phi[i][j2] = phi[i][j2 + 1];
    }
}
/***** End bcforphi *****/


/*****  poiseq : solve poisson's equation of pressure field *****/
void SolveFlow::poiseq(double &dt, vector2double &u, vector2double &v, vector2double &p) {
    // compute (RHS)
    vector2double rhs(u.size(), std::vector<double>(u.front().size()));

    for (int i = 1; i < mx - 1; i++) {
        for (int j = 1; j < my - 1; j++) {
            if ((i1 < i && i < i2) && (j1 < j && j < j2)) {
                continue;
            }
            double ux = 0.5 * dxi * (u[i + 1][j] - u[i - 1][j]);  // dxi=1/dx
            double uy = 0.5 * dyi * (u[i][j + 1] - u[i][j - 1]);  // dyi=1/dy
            double vx = 0.5 * dxi * (v[i + 1][j] - v[i - 1][j]);  // dxi=1/dx
            double vy = 0.5 * dyi * (v[i][j + 1] - v[i][j - 1]);  // dyi=1/dy
            rhs[i][j] = 1.0 / dt * (ux + vy) - (ux * ux + 2.0 * uy * vx + vy * vy);
        }
    }

    // iterations
    for (int itr = 0; itr <= maxitp; itr++) {  // iterations
        double res = 0.0;


        for (int i = 1; i < mx - 1; i++) {
            for (int j = 1; j < my - 1; j++) {
                // ignore inside rectangle
                if ((i1 < i && i < i2) && (j1 < j && j < j2)) {
                    continue;
                }
                double dp = dxi2 * (p[i + 1][j] + p[i - 1][j]) +
                            dyi2 * (p[i][j + 1] + p[i][j - 1]) - rhs[i][j];
                dp = dp * dxdy - p[i][j];
                res += dp * dp;          // residual
                p[i][j] += omegap * dp;  // when omegap > 1, accerelate convergence
            }
        }

        bcforp(p);  // set BC
        res = std::sqrt(1.0 / double(mx * my) * res);
        //    cout<<itr<<" "<<res<<endl;

        if (res < errorp || itr == maxitp) {
            // int resp = res;
            // int itrp = itr;
            // std::cout << resp << " " << itrp << std::endl;
            break;
        }
    }
    //    cout<<itrp<<" "<<resp<<endl;
}

/***** veloeq: solve velocity field: kawamura scheme *****/
void SolveFlow::veloeq(double &dt, vector2double &u, vector2double &v, vector2double &p) {
    vector2double urhs(u.size(), std::vector<double>(u.front().size())),
        vrhs(v.size(), std::vector<double>(v.front().size()));
    // pressure gradient
    for (int i = 1; i < mx - 1; i++) {
        for (int j = 1; j < my - 1; j++) {
            if ((i1 < i && i < i2) && (j1 < j && j < j2)) {
                continue;
            }
            urhs[i][j] = -0.5 * dxi * (p[i + 1][j] - p[i - 1][j]);
            vrhs[i][j] = -0.5 * dyi * (p[i][j + 1] - p[i][j - 1]);
        }
    }

    // Viscous Term
    for (int i = 1; i < mx - 1; i++) {
        for (int j = 1; j < my - 1; j++) {
            if ((i1 < i && i < i2) && (j1 < j && j < j2)) {
                continue;
            }
            urhs[i][j] += rei * (dxi2 * (u[i + 1][j] - 2.0 * u[i][j] + u[i - 1][j]) +
                                 dyi2 * (u[i][j + 1] - 2.0 * u[i][j] + u[i][j - 1]));
            vrhs[i][j] += rei * (dxi2 * (v[i + 1][j] - 2.0 * v[i][j] + v[i - 1][j]) +
                                 dyi2 * (v[i][j + 1] - 2.0 * v[i][j] + v[i][j - 1]));
        }
    }

    // Advection term in x
    // calculate points needed for kuwahara scheme
    for (int j = j1 + 1; j <= j2 - 1; j++) {
        u[i1 + 1][j] = 2.0 * u[i1][j] - u[i1 - 1][j];
        u[i2 - 1][j] = 2.0 * u[i2][j] - u[i2 + 1][j];
        v[i1 + 1][j] = 2.0 * v[i1][j] - v[i1 - 1][j];
        v[i2 - 1][j] = 2.0 * v[i2][j] - v[i2 + 1][j];
    }

    for (int i = 2; i < mx - 2; i++) {
        for (int j = 0; j < my; j++) {
            if ((i1 < i && i < i2) && (j1 < j && j < j2)) {  // ignore inside rectangle
                continue;
            }
            urhs[i][j] += -1.0 / 12 * u[i][j] * dxi *
                              (-u[i + 2][j] + 8.0 * (u[i + 1][j] - u[i - 1][j]) + u[i - 2][j]) -
                          0.25 * dxi * abs(u[i][j]) *
                              (u[i + 2][j] - 4.0 * u[i + 1][j] + 6.0 * u[i][j] - 4.0 * u[i - 1][j] +
                               u[i - 2][j]);
            vrhs[i][j] += -1.0 / 12 * u[i][j] * dxi *
                              (-v[i + 2][j] + 8.0 * (v[i + 1][j] - v[i - 1][j]) + v[i - 2][j]) -
                          0.25 * dxi * abs(u[i][j]) *
                              (v[i + 2][j] - 4.0 * v[i + 1][j] + 6.0 * v[i][j] - 4.0 * v[i - 1][j] +
                               v[i - 2][j]);
        }
    }

    // advection term in y
    for (int i = i1 + 1; i <= i2 - 1; i++) {
        u[i][j1 + 1] = 2.0 * u[i][j1] - u[i][j1 - 1];
        u[i][j2 - 1] = 2.0 * u[i][j2] - u[i][j2 + 1];
        v[i][j1 + 1] = 2.0 * v[i][j1] - v[i][j1 - 1];
        v[i][j2 - 1] = 2.0 * v[i][j2] - v[i][j2 + 1];
    }

    for (int i = 0; i < mx; i++) {
        for (int j = 2; j < my - 2; j++) {
            if ((i1 < i && i < i2) && (j1 < j && j < j2)) {
                continue;
            }
            urhs[i][j] += -1.0 / 12 * v[i][j] * dyi *
                              (-u[i][j + 2] + 8.0 * (u[i][j + 1] - u[i][j - 1]) + u[i][j - 2]) -
                          0.25 * dyi * abs(v[i][j]) *
                              (u[i][j + 2] - 4.0 * u[i][j + 1] + 6.0 * u[i][j] - 4.0 * u[i][j - 1] +
                               u[i][j - 2]);
            vrhs[i][j] += -1.0 / 12 * v[i][j] * dyi *
                              (-v[i][j + 2] + 8.0 * (v[i][j + 1] - v[i][j - 1]) + v[i][j - 2]) -
                          0.25 * dyi * abs(v[i][j]) *
                              (v[i][j + 2] - 4.0 * v[i][j + 1] + 6.0 * v[i][j] - 4.0 * v[i][j - 1] +
                               v[i][j - 2]);
        }
    }

    // update
    for (int i = 1; i < mx - 1; i++) {
        for (int j = 1; j < my - 1; j++) {
            if ((i1 < i && i < i2) && (j1 < j && j < j2)) {
                continue;
            }
            u[i][j] += dt * urhs[i][j];
            v[i][j] += dt * vrhs[i][j];
        }
    }
}

/***** veloeq: solve velocity field: kawamura scheme *****/
void SolveFlow::denseq(double &dt, vector2double &u, vector2double &v, vector2double &phi) {
    vector2double phirhs(phi.size(), std::vector<double>(phi.front().size()));
    // Diffusion Term
    for (int i = 1; i < mx - 1; i++) {
        for (int j = 1; j < my - 1; j++) {
            if ((i1 < i && i < i2) && (j1 < j && j < j2)) {
                continue;  // skip inside prism
            }
            phirhs[i][j] += d * (dxi2 * (phi[i + 1][j] - 2.0 * phi[i][j] + phi[i - 1][j]) +
                                 dyi2 * (phi[i][j + 1] - 2.0 * phi[i][j] + phi[i][j - 1]));
        }
    }

    // Advection term in x
    // calculate points needed for kuwahara scheme
    for (int j = j1 + 1; j <= j2 - 1; j++) {
        phi[i1 + 1][j] = 2.0 * phi[i1][j] - phi[i1 - 1][j];
        phi[i2 - 1][j] = 2.0 * phi[i2][j] - phi[i2 + 1][j];
    }

    for (int i = 2; i < mx - 2; i++) {
        for (int j = 0; j < my; j++) {
            if ((i1 < i && i < i2) && (j1 < j && j < j2)) {  // ignore inside rectangle
                continue;
            }
            phirhs[i][j] +=
                -1.0 / 12 * u[i][j] * dxi *
                    (-phi[i + 2][j] + 8.0 * (phi[i + 1][j] - phi[i - 1][j]) + phi[i - 2][j]) -
                0.25 * dxi * abs(u[i][j]) *
                    (phi[i + 2][j] - 4.0 * phi[i + 1][j] + 6.0 * phi[i][j] - 4.0 * phi[i - 1][j] +
                     phi[i - 2][j]);
        }
    }

    // advection term in y
    for (int i = i1 + 1; i <= i2 - 1; i++) {
        phi[i][j1 + 1] = 2.0 * phi[i][j1] - phi[i][j1 - 1];
        phi[i][j2 - 1] = 2.0 * phi[i][j2] - phi[i][j2 + 1];
    }

    for (int i = 0; i < mx; i++) {
        for (int j = 2; j < my - 2; j++) {
            if ((i1 < i && i < i2) && (j1 < j && j < j2)) {
                continue;
            }
            phirhs[i][j] +=
                -1.0 / 12 * v[i][j] * dyi *
                    (-phi[i][j + 2] + 8.0 * (phi[i][j + 1] - phi[i][j - 1]) + phi[i][j - 2]) -
                0.25 * dyi * abs(v[i][j]) *
                    (phi[i][j + 2] - 4.0 * phi[i][j + 1] + 6.0 * phi[i][j] - 4.0 * phi[i][j - 1] +
                     phi[i][j - 2]);
        }
    }

    // update
    for (int i = 1; i < mx - 1; i++) {
        for (int j = 1; j < my - 1; j++) {
            if ((i1 < i && i < i2) && (j1 < j && j < j2)) {
                continue;
            }
            phi[i][j] += dt * phirhs[i][j];
        }
    }
}
