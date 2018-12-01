CXX = g++
TARGET = boundarylayer
CXXFLAGS = -Wall -O2 -std=c++11
SRCS = ./src/boundary_layer.cpp
OBJS := $(SRCS:.cpp=.o) #SRCSの各ファイルのサフィックスの.cppを.oに変換

$(TARGET): $(OBJS)
	$(CXX) -o $@ $(OBJS)

clean:
	rm -f $(OBJS)
