#Makefile for c++11

#プログラム名とオブジェクトファイル名
TARGET = ic_ns2d
OBJS = ./src/ic_ns2d.cpp ./src/setgrid.cpp ./src/solveflow.cpp ./src/streamline.cpp

#定義済みマクロの再定義
CXX = g++
CXXFLAGS = -Wall -g -ggdb -O0 -std=c++11

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
ic_ns2d.o:  ./src/flowparam.hpp ./src/setgrid.hpp ./src/streamline.hpp ./src/solveflow.hpp
setgrid.o : ./src/setgrid.hpp
solveflow.o : ./src/flowparam.hpp ./src/setgrid.hpp ./src/solveflow.hpp
streamline.o: ./src/flowparam.hpp ./src/setgrid.hpp ./src/streamline.hpp
