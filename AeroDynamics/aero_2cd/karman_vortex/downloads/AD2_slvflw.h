#ifndef AD2_SlvFlw_H
#define AD2_SlvFlw_H

using namespace std;

const double re    = 70;
double       rei   = 1 / re;
int          nlast = 5000, nlp = 50;

class SlvFlw {
public:
    static void slvflw(double& dt, vector<vector<double>>& x, vector<vector<double>>& y,
                       vector<vector<double>>& u, vector<vector<double>>& v,
                       vector<vector<double>>& p);
    static void intcnd(double& nbegin, double& time, vector<vector<double>>& u,
                       vector<vector<double>>& v, vector<vector<double>>& p);
    static void bcforp(vector<vector<double>>& p);
    static void bcforv(vector<vector<double>>& u, vector<vector<double>>& v);
    static void poiseq(double& dt, vector<vector<double>>& u, vector<vector<double>>& v,
                       vector<vector<double>>& p);
    static void veloeq(double& dt, vector<vector<double>>& u, vector<vector<double>>& v,
                       vector<vector<double>>& p);

private:
    static const int mx = SetGrd::mx, my = SetGrd::my, i1 = SetGrd::i1, i2 = SetGrd::i2,
                     j1 = SetGrd::j1, j2 = SetGrd::j2;  // copy values, dxi:1/dx
    static constexpr double dxi = 1.0 / SetGrd::dx, dyi = 1.0 / SetGrd::dy, dxi2 = dxi * dxi,
                            dyi2 = dyi * dyi,
                            dxdy = 1.0 / (dxi * dxi + dyi * dyi) * 0.5,  //=2/(1/dx^2+1/dy^2)
        omegap = 1.1, errorp = 1e-4;
    static const int maxitp = 200;
    static int       itrp;
    static double    resp;
};
int    SlvFlw::itrp = 0;
double SlvFlw::resp = 1.0;


/***** slvflw  *****/
void SlvFlw::slvflw(double& dt, vector<vector<double>>& x, vector<vector<double>>& y,
                    vector<vector<double>>& u, vector<vector<double>>& v,
                    vector<vector<double>>& p) {
    fstream fout;
    fout.open("pmap.txt", ios::out);

    double nbegin, time;
    cout << " slvflw" << endl;
    intcnd(nbegin, time, u, v, p);
    bcforp(p);
    bcforv(u, v);

    //  cout<<" "<<dxi<<" "<<dxi2<<endl;
    cout << "Step/resp/itrp/cd/cl/cp1/cp2" << endl;


    for (int n = 1; n <= nlast; n++) {
        double nstep = n + nbegin;
        time += dt;

        poiseq(dt, u, v, p);
        bcforp(p);

        veloeq(dt, u, v, p);
        bcforv(u, v);

        double cd = 0.0;
        for (int j = j1; j <= j2 - 1; j++) {
            double cpfore = 1.0 * (p[i1][j] + p[i1][j + 1]);
            double cpback = 1.0 * (p[i2][j] + p[i2][j + 1]);
            cd += SetGrd::dy * (cpfore - cpback);
        }

        double cl = 0.0;
        for (int i = i1; i <= i2 - 1; i++) {
            double cpbtm = 1.0 * (p[i][j1] + p[i + 1][j1]);
            double cptop = 1.0 * (p[i][j2] + p[i + 1][j2]);
            cl += SetGrd::dx * (cpbtm - cptop);
        }

        if (!(n % nlp)) {
            double cp1 = 2.0 * p[i2 + i2 - i1][j1];
            double cp2 = 2.0 * p[i2 + i2 - i1][j2];
            cout << nstep << " " << resp << " " << itrp << " " << cd << " " << cl << " " << cp1
                 << " " << cp2 << endl;
        }

        if (n == 3500) {
            for (int i = 1; i <= mx; i++) {
                for (int j = 1; j <= my; j++) {
                    fout << x[i][j] << " " << y[i][j] << " " << 2.0 * p[i][j] << endl;
                }
                fout << endl;
            }
        }
    }

    fout.close();
}
/***** End slvflw  *****/

/***** intcnd  *****/
void SlvFlw::intcnd(double& nbegin, double& time, vector<vector<double>>& u,
                    vector<vector<double>>& v, vector<vector<double>>& p) {
    nbegin = 0;
    time   = 0.0;

    for (int i = 1; i <= mx; i++) {
        for (int j = 1; j <= my; j++) {
            u[i][j] = 1.0;
            v[i][j] = 0.0;
            p[i][j] = 0.0;
        }
    }

    cout << " intcnd" << endl;
}
/***** End intcnd  *****/

/***** bcforp  *****/
void SlvFlw::bcforp(vector<vector<double>>& p) {
    //  cout<<" bcforphoge";
    for (int j = 1; j <= my; j++) {
        int i   = 1;
        p[i][j] = 0.0;  // inflow

        i       = mx;
        p[i][j] = 0.0;  // downstream
    }

    for (int i = 1; i <= mx; i++) {
        int j   = 1;
        p[i][j] = 0.0;  // bottom

        j       = my;
        p[i][j] = 0.0;  // top
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
/***** End bcforp  *****/


/***** bcforv  *****/
void SlvFlw::bcforv(vector<vector<double>>& u, vector<vector<double>>& v) {
    for (int j = 1; j <= my; j++) {
        int i       = 1;
        u[i][j]     = 1.0;  // inflow
        v[i][j]     = 0.0;
        u[i - 1][j] = 1.0;
        v[i - 1][j] = 0.0;

        i           = mx;
        u[i][j]     = 2.0 * u[i - 1][j] - u[i - 2][j];  // downstream
        v[i][j]     = 2.0 * v[i - 1][j] - v[i - 2][j];
        u[i + 1][j] = 2.0 * u[i][j] - u[i - 1][j];
        v[i + 1][j] = 2.0 * v[i][j] - v[i - 1][j];
    }

    for (int i = 1; i <= mx; i++) {
        int j       = 1;
        u[i][j]     = 2.0 * u[i][j + 1] - u[i][j + 2];  // bottom
        v[i][j]     = 2.0 * v[i][j + 1] - v[i][j + 2];
        u[i][j - 1] = 2.0 * u[i][j] - u[i][j + 1];
        v[i][j - 1] = 2.0 * v[i][j] - v[i][j + 1];

        j           = my;
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


/***** poiseq *****/
void SlvFlw::poiseq(double& dt, vector<vector<double>>& u, vector<vector<double>>& v,
                    vector<vector<double>>& p) {
    vector<vector<double>> rhs(u.size(), vector<double>(u.front().size()));

    for (int i = 2; i <= mx - 1; i++) {  // compute RHS
        for (int j = 2; j <= my - 1; j++) {
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

    for (int itr = 0; itr <= maxitp; itr++) {  // iterations
        double res = 0.0;

        for (int i = 2; i <= mx - 1; i++) {
            for (int j = 2; j <= my - 1; j++) {
                if ((i1 < i && i < i2) && (j1 < j && j < j2)) {
                    continue;
                }
                double dp = dxi2 * (p[i + 1][j] + p[i - 1][j]) +
                            dyi2 * (p[i][j + 1] + p[i][j - 1]) - rhs[i][j];
                dp = dp * dxdy - p[i][j];
                res += dp * dp;
                p[i][j] += omegap * dp;
            }
        }

        bcforp(p);  // set BC
        res = sqrt(1.0 / double(mx * my) * res);
        //    cout<<itr<<" "<<res<<endl;

        if (res < errorp || itr == maxitp) {
            resp = res;
            itrp = itr;
            break;
        }
    }
    //    cout<<itrp<<" "<<resp<<endl;
}
/***** End poiseq *****/


/***** veloeq *****/
void SlvFlw::veloeq(double& dt, vector<vector<double>>& u, vector<vector<double>>& v,
                    vector<vector<double>>& p) {
    vector<vector<double>> urhs(u.size(), vector<double>(u.front().size())),
        vrhs(v.size(), vector<double>(v.front().size()));
    for (int i = 2; i <= mx - 1; i++) {  // pressure gradient
        for (int j = 2; j <= my - 1; j++) {
            if ((i1 < i && i < i2) && (j1 < j && j < j2)) {
                continue;
            }
            urhs[i][j] = -0.5 * dxi * (p[i + 1][j] - p[i - 1][j]);
            vrhs[i][j] = -0.5 * dyi * (p[i][j + 1] - p[i][j - 1]);
        }
    }

    for (int i = 2; i <= mx - 1; i++) {  // Viscous Term
        for (int j = 2; j <= my - 1; j++) {
            if ((i1 < i && i < i2) && (j1 < j && j < j2)) {
                continue;
            }
            urhs[i][j] += rei * (dxi2 * (u[i + 1][j] - 2.0 * u[i][j] + u[i - 1][j]) +
                                 dyi2 * (u[i][j + 1] - 2.0 * u[i][j] + u[i][j - 1]));
            vrhs[i][j] += rei * (dxi2 * (v[i + 1][j] - 2.0 * v[i][j] + v[i - 1][j]) +
                                 dyi2 * (v[i][j + 1] - 2.0 * v[i][j] + v[i][j - 1]));
        }
    }

    for (int j = j1 + 1; j <= j2 - 1; j++) {  // advection term in x
        u[i1 + 1][j] = 2.0 * u[i1][j] - u[i1 - 1][j];
        u[i2 - 1][j] = 2.0 * u[i2][j] - u[i2 + 1][j];
        v[i1 + 1][j] = 2.0 * v[i1][j] - v[i1 - 1][j];
        v[i2 - 1][j] = 2.0 * v[i2][j] - v[i2 + 1][j];
    }

    for (int i = 2; i <= mx - 1; i++) {
        for (int j = 2; j <= my - 1; j++) {
            if ((i1 < i && i < i2) && (j1 < j && j < j2)) {
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

    for (int i = i1 + 1; i <= i2 - 1; i++) {  // advection term in y
        u[i][j1 + 1] = 2.0 * u[i][j1] - u[i][j1 - 1];
        u[i][j2 - 1] = 2.0 * u[i][j2] - u[i][j2 + 1];
        v[i][j1 + 1] = 2.0 * v[i][j1] - v[i][j1 - 1];
        v[i][j2 - 1] = 2.0 * v[i][j2] - v[i][j2 + 1];
    }

    for (int i = 2; i <= mx - 1; i++) {
        for (int j = 2; j <= my - 1; j++) {
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

    for (int i = 2; i <= mx - 1; i++) {  // update
        for (int j = 2; j <= my - 1; j++) {
            if ((i1 < i && i < i2) && (j1 < j && j < j2)) {
                continue;
            }
            u[i][j] += dt * urhs[i][j];
            v[i][j] += dt * vrhs[i][j];
        }
    }
}
/***** End veloeq *****/

#endif
