#Makefile for c++11

#プログラム名とオブジェクトファイル名
TARGET = report2
OBJS = ../src/report2.cpp ../src/non_linear.cpp ../src/kalman_filter.cpp

#定義済みマクロの再定義
CXX = g++
CXXFLAGS = -Wall -g -ggdb -std=c++11

#サフィックスルール適用対象の拡張子の定義
.SUFFIXES: .cpp .o  #SRCSの各ファイルのサフィックスの.cppを.oに変換

#プライマリーターゲット
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o  $(TARGET) $^

#サフィックスルール
.cpp.o:
	$(CXX) $(CXXFLAGS) $@ -c $<

#ファイル削除用ターゲット
.PHONY: clean
clean:
	rm -f $(OBJS)

#ヘッダーファイルの依存関係
report2.o:  ../src/non_linear.hpp ../src/kalman_filter.hpp \
	 ../src/Eigen/Core ../src/Eigen/Geometry ../src/unsupported/Eigen/MatrixFunctions
non_linear.o: ../src/non_linear.hpp ../src/Eigen/Core ../src/Eigen/Geometry
kalman_filter.o: ../src/kalman_filter.hpp ../src/Eigen/Core  \
	../src/Eigen/Geometry ../src/unsupported/Eigen/MatrixFunctions ../src/non_linear.hpp
