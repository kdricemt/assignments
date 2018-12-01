#ifndef FLOWPARAM_H
#define FLOWPARAM_H

struct FlowParam {
    static constexpr double re    = 70.0;  // reynolds number
    static constexpr double rei   = 1 / re;
    static constexpr double cfl   = 0.2;     // cfl number
    static constexpr int    nlast = 5000;    // number of time steps
    static constexpr int    nlp   = 10;      // interval of steps for calculation
    static constexpr double d     = 1.0e-5;  //拡散係数
};

#endif
