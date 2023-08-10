
#
# Please Read - This makefile for the parallel RLA implementation
#


CXX		 = g++
CXXFLAGS = -Wall -std=c++17 -g -pg
TARGET	 = RLA_parallel
SRC		 = RLA_parallel.cpp
LIB      = -lbenchmark -lpthread
$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(SRC) $(LIB) -o $(TARGET)

.PHONY: clean
clean:
	rm -f $(TARGET)

