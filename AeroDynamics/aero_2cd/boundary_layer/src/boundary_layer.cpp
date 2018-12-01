#include <chrono>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <vector>

class BoundaryLayer {
private:
    //気体パラメータ
    const double rou0_  = 1.225;                                           //密度
    const double t0_    = 293.5;                                           //温度
    const double visc_  = 1.458e-06 * std::pow(t0_, 1.5) / (t0_ + 110.4);  //粘性係数
    const double width_ = 0.5;                                             //板の長さ
    const double v0_    = 20.0;                                            //流速
    const double re0_   = rou0_ * v0_ * width_ / visc_;                    //レイノルズ数
    const double blt_   = width_ * 5.3 / std::sqrt(re0_);                  //境界層厚さ(m)

    //格子パラメータ
    const int           nx_   = 50000;  // x軸方向格子数
    const double        dx_   = width_ / nx_;
    const int           ny_   = 200;  // y軸方向格子数
    const double        ymax_ = blt_ * 2.0;
    const double        dy_   = ymax_ / ny_;
    std::vector<double> vx_;  // x座標
    std::vector<double> vy_;  // y座標

    //初期条件
    const double atm_ = 101325;
    const double vwall_;  //壁での吸い込み／湧き出し
    const double dpdx_;   //圧力勾配

    //流速
    std::vector<std::vector<double>> vu_;  // x軸方向速さ
    std::vector<std::vector<double>> vv_;  // y軸方向速さ

    //境界層厚さ
    std::vector<double> bt_;  //境界層厚さ
    std::vector<double> dt_;  //排除厚さ
    std::vector<double> mt_;  //運動量厚さ
    std::vector<double> et_;  //エネルギー厚さ

    //摩擦係数
    std::vector<double> cf_;

    // Shape factor
    std::vector<double> h12_;  //形状係数 H12 = bt_/mt_
    std::vector<double> h23_;  //形状係数 H23 = mt_/et_

public:
    BoundaryLayer(double vwall, double dpdx_atm);
    int cal_grids();
    int write_files();
    int plot_boundary_layer();
    int plot_shape_factor();
};

BoundaryLayer::BoundaryLayer(double vwall, double dpdx_atm)
    : vwall_(vwall), dpdx_(dpdx_atm * atm_) {
    std::printf("========================================================\n");
    std::printf("BOUNDARY LAYER CFD\n");
    std::printf("========================================================\n");
    std::printf("  Number of grids(x,y)             :    %d x %d\n", nx_, ny_);
    std::printf("  Air Speed outside BoundaryLayer  :    %.1f [m/s]\n", v0_);
    std::printf("  Divergence at wall               :    %.3f [m/s]\n", vwall_);
    std::printf("  Pressure gradient                :    %.1f [Pa/m]\n", dpdx_);
    std::printf("  Reynolds number                  :    %.1f \n", re0_);
    std::printf("  Boundarylayer thickness(Blaudius):    %.1f [mm]\n", blt_ * 1000);
    std::printf("========================================================\n");


    std::vector<double> vu_1;
    std::vector<double> vv_1;
    //格子条件,スタート条件
    vx_.reserve(nx_);
    vy_.reserve(ny_);
    vu_1.reserve(ny_);
    vv_1.reserve(ny_);
    for (int column = 0; column < nx_; column++) {
        vx_.push_back(dx_ * column);
    }
    for (int row = 0; row < ny_; row++) {
        vy_.push_back(dy_ * row);
        // 1列目の条件を入力
        vu_1.push_back(v0_);
        vv_1.push_back(0.0);
    }
    vu_.push_back(vu_1);
    vv_.push_back(vv_1);
    bt_.push_back(0.0);
    dt_.push_back(0.0);
    mt_.push_back(0.0);
    et_.push_back(0.0);
    h12_.push_back(0.0);
    h23_.push_back(0.0);
    // x=0とすると無限大に発散するので、適当に小さい値(dx_/5)で置き換える
    cf_.push_back(0.664 / std::sqrt(rou0_ * v0_ * (dx_ / 5) / visc_));
}

//前の列(x軸方向)の結果を利用して、次の列のvu,vv(vu_new,vv_new)を計算する
int BoundaryLayer::cal_grids() {
    auto start = std::chrono::system_clock::now();  // 計測スタート時刻を保存
    // 高速化のために、あらかじめメモリを確保しておく
    vu_.reserve(nx_ * ny_);
    vv_.reserve(nx_ * ny_);

    for (int x_column = 1; x_column < nx_; x_column++) {
        std::vector<double> vu_new;
        std::vector<double> vv_new;
        vu_new.reserve(ny_);
        vv_new.reserve(ny_);
        double u_e = v0_;  //境界層外部の流速

        //この列のvu_を計算
        vu_new.push_back(0.0);
        for (int y_row = 1; y_row < ny_ - 1; y_row++) {
            double dudy1 =
                (vu_[x_column - 1][y_row + 1] - vu_[x_column - 1][y_row - 1]) / (2.0 * dy_);
            double dudy2 = (vu_[x_column - 1][y_row + 1] - 2.0 * vu_[x_column - 1][y_row] +
                            vu_[x_column - 1][y_row - 1]) /
                           (dy_ * dy_);
            double dudy3 = -(1 / rou0_) * dpdx_;
            double dudx  = (dudy2 * visc_ / rou0_ - vv_[x_column - 1][y_row] * dudy1 + dudy3) /
                          vu_[x_column - 1][y_row];
            vu_new.push_back(vu_[x_column - 1][y_row] + dx_ * dudx);
        }

        if (dpdx_ != 0) u_e = vu_new[ny_ - 2];  //圧力勾配がある場合はu_e=v0_ではない
        vu_new.push_back(u_e);

        //この列のvvを計算
        vv_new.push_back(vwall_);
        for (int y_row = 1; y_row < ny_ - 1; y_row++) {
            double dudx1 = (vu_new[y_row] - vu_[x_column - 1][y_row]) / dx_;
            double dudx2 = (vu_new[y_row + 1] - vu_[x_column - 1][y_row + 1]) / dx_;
            vv_new.push_back(vv_new[y_row - 1] - (dudx1 + dudx2) * dy_ / 2.0);
        }

        //この列のbtを計算
        double velocity_criteria = 0.995 * u_e;
        double y_index           = 0;
        double bt_new            = 0;
        while (vu_new[y_index] < velocity_criteria) {
            y_index++;
        }
        //線形近似
        bt_new = (vy_[y_index - 1] * (vu_new[y_index] - velocity_criteria) +
                  vy_[y_index] * (velocity_criteria - vu_new[y_index - 1])) /
                 (vu_new[y_index] - vu_new[y_index - 1]);

        bt_.push_back(bt_new);

        //この列のdtを計算
        double dt_new = 0;
        for (int y_row = 0; y_row < ny_; y_row++) {
            dt_new += (1 - vu_new[y_row] / u_e) * dy_;
        }
        dt_.push_back(dt_new);

        //この列のmtを計算
        double mt_new = 0;
        for (int y_row = 0; y_row < ny_; y_row++) {
            mt_new += (vu_new[y_row] / u_e) * (1 - vu_new[y_row] / u_e) * dy_;
        }
        mt_.push_back(mt_new);

        //この列のetを計算
        double et_new = 0;
        for (int y_row = 0; y_row < ny_; y_row++) {
            et_new += (vu_new[y_row] / u_e) * (1 - std::pow(vu_new[y_row] / u_e, 2)) * dy_;
        }
        et_.push_back(et_new);

        //この列のh12,h23を計算
        h12_.push_back(dt_new / mt_new);
        h23_.push_back(mt_new / et_new);

        //この列のCfを計算
        double tau_new = 0;
        tau_new        = visc_ * (vu_new[1] - vu_new[0]) / dy_;
        cf_.push_back(2 * tau_new / (rou0_ * u_e * u_e));

        //結合
        vu_.push_back(vu_new);
        vv_.push_back(vv_new);

        // vu_new,vv_newのメモリを解放
        vu_new.clear();
        vv_new.clear();
        vu_new.shrink_to_fit();
        vv_new.shrink_to_fit();
    }

    auto end  = std::chrono::system_clock::now();  // 計測終了時刻を保存
    auto dur  = end - start;                       // 要した時間を計算
    auto msec = std::chrono::duration_cast<std::chrono::milliseconds>(dur).count();
    // 要した時間をミリ秒（1/1000秒）に変換して表示
    std::cout << "  calulation time                  :     " << msec << " [msec]\n";
    std::printf("  BoundaryLayer Thickness at edge  :     %.2f [mm]\n", bt_[nx_ - 1] * 1000);
    std::printf("  Displacement Thickness at edge   :     %.2f [mm]\n", dt_[nx_ - 1] * 1000);
    std::printf("  Momentum Thickness at edge       :     %.2f [mm]\n", mt_[nx_ - 1] * 1000);
    std::printf("  Energy Thickness at edge         :     %.2f [mm]\n", et_[nx_ - 1] * 1000);
    std::printf("  H12 at edge                      :     %.2f \n", h12_[nx_ - 1]);
    std::printf("  H23 at edge                      :     %.2f \n", h23_[nx_ - 1]);
    return 0;
}

int BoundaryLayer::write_files() {
    //塗りつぶし用速度vectorデータ書き込み
    FILE *      fid_flow_vector_all;
    const char *data_ffva = {"./data_files/flow_vector_all.txt"};
    fid_flow_vector_all   = fopen(data_ffva, "w");
    int nx_devide         = 10;  // vectorを何列表示するか
    for (int x = 0; x < nx_; x += int(nx_ / nx_devide)) {
        for (int y = 0; y < ny_; y++) {
            fprintf(fid_flow_vector_all, "%f\t%f\t%f\t%f\n", vx_[x], vy_[y], vu_[x][y], vv_[x][y]);
        }
        fprintf(fid_flow_vector_all, "\n");
    }

    //流���用��度vectorデータ書き込み
    FILE *      fid_flow_vector_part;
    const char *data_ffvp = {"./data_files/flow_vector_part.txt"};
    fid_flow_vector_part  = fopen(data_ffvp, "w");
    int ny_devide         = 20;  // vectorを一列ごとに何行表示するか
    for (int x = 0; x < nx_; x += int(nx_ / nx_devide)) {
        for (int y = 0; y < ny_; y += int(ny_ / ny_devide)) {
            fprintf(fid_flow_vector_part, "%f\t%f\t%f\t%f\n", vx_[x], vy_[y], vu_[x][y], vv_[x][y]);
        }
        fprintf(fid_flow_vector_part, "\n");
    }

    //境界層厚�����データの書き込み
    FILE *      fid_bt;
    const char *data_bt = {"./data_files/bt.txt"};
    fid_bt              = fopen(data_bt, "w");
    for (int x = 0; x < nx_; x++) {
        fprintf(fid_bt, "%f\t%f\n", vx_[x], bt_[x]);
    }

    //排除厚さデータの書き込み
    FILE *      fid_dt;
    const char *data_dt = {"./data_files/dt.txt"};
    fid_dt              = fopen(data_dt, "w");
    for (int x = 0; x < nx_; x++) {
        fprintf(fid_dt, "%f\t%f\n", vx_[x], dt_[x]);
    }

    //運動量厚さデータの書き込み
    FILE *      fid_mt;
    const char *data_mt = {"./data_files/mt.txt"};
    fid_mt              = fopen(data_mt, "w");
    for (int x = 0; x < nx_; x++) {
        fprintf(fid_mt, "%f\t%f\n", vx_[x], mt_[x]);
    }

    //エネルギー厚さデータの書き込み
    FILE *      fid_et;
    const char *data_et = {"./data_files/et.txt"};
    fid_et              = fopen(data_et, "w");
    for (int x = 0; x < nx_; x++) {
        fprintf(fid_et, "%f\t%f\n", vx_[x], et_[x]);
    }

    // Cfデータの書き込み
    FILE *      fid_cf;
    const char *data_cf = {"./data_files/cf.txt"};
    fid_cf              = fopen(data_cf, "w");
    for (int x = 0; x < nx_; x++) {
        fprintf(fid_cf, "%f\t%f\n", vx_[x], cf_[x]);
    }

    // H12,H23データの書き込み 条件により書き込むファイルを変える
    FILE *      fid_h12;
    FILE *      fid_h23;
    const char *data_h12;
    const char *data_h23;
    if (vwall_ > 0) {
        data_h12 = {"./data_files/h12_div_plus.txt"};
        data_h23 = {"./data_files/h23_div_plus.txt"};
    } else if (vwall_ < 0) {
        data_h12 = {"./data_files/h12_div_minus.txt"};
        data_h23 = {"./data_files/h23_div_minus.txt"};
    } else if (dpdx_ > 0) {
        data_h12 = {"./data_files/h12_dpdx_plus.txt"};
        data_h23 = {"./data_files/h23_dpdx_plus.txt"};
    } else if (dpdx_ < 0) {
        data_h12 = {"./data_files/h12_dpdx_minus.txt"};
        data_h23 = {"./data_files/h23_dpdx_minus.txt"};
    } else {
        data_h12 = {"./data_files/h12_standard.txt"};
        data_h23 = {"./data_files/h23_standard.txt"};
    }
    fid_h12 = fopen(data_h12, "w");
    fid_h23 = fopen(data_h23, "w");
    for (int x = 0; x < nx_; x++) {
        fprintf(fid_h12, "%f\t%f\n", vx_[x], h12_[x]);
        fprintf(fid_h23, "%f\t%f\n", vx_[x], h23_[x]);
    }
    return 0;
}

int BoundaryLayer::plot_boundary_layer() {
    FILE *gid = popen("gnuplot", "w");
    fprintf(gid, "load 'boundary_layer.gp'");
    pclose(gid);
    return 0;
}

int BoundaryLayer::plot_shape_factor() {
    FILE *gid = popen("gnuplot", "w");
    fprintf(gid, "load 'shape_data.gp'");
    pclose(gid);
    return 0;
}

int main() {
    BoundaryLayer b1(-0.015, -0.00);
    b1.cal_grids();
    b1.write_files();
    // b1.plot_boundary_layer();
    // b1.plot_shape_factor();
    return 0;
};
