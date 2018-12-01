#ifndef SETGRID_H
#define SETGRID_H

#include <vector>

typedef std::vector<std::vector<double>> vector2double;

// Set Grid Parameters
/* ひとつはROM化可能な値、もうひとつは実行時にしか決まらないがいったん初期化したあとは二度と変更されない値である。
C++11以降、前者はconstexprが受け持ち、後者はconstが受け持つことになった。*/
class SetGrid {
public:
    // i1,i2,j1,j2 : wall edge
    static constexpr int    mx = 400, i1 = 95, i2 = 105, my = 200, j1 = 95, j2 = 105;
    static constexpr double dx = 1 / double(i2 - i1), dy = 1 / double(j2 - j1);
    static double           setgrd(vector2double &x, vector2double &y);
};

#endif
