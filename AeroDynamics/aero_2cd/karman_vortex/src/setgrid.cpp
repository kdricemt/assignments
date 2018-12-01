#include <cstdio>
#include <vector>
#include "setgrid.hpp"
// setgrid 計算格子の設定
using std::vector;

double SetGrid::setgrd(vector<vector<double>>& x, vector<vector<double>>& y) {
    double icent = 0.5 * (i1 + i2), jcent = 0.5 * (j1 + j2);
    for (int i = 0; i < mx; i++) {
        for (int j = 0; j < my; j++) {
            x[i][j] = dx * double(i - icent);
            y[i][j] = dy * double(j - jcent);
        }
    }

    if (dx >= dy) {
        return dx;
    } else {
        return dy;
    }
}
