CC=gcc
CXX=g++
CFLAGS=-O2
CXXFLAGS=-O2 --std=c++11
LFLAGS=-s -lgdi32

all: viewbmp mandelbrot nqueens sierpinski tree cantor hilbert

viewbmp: viewbmp.cpp
	$(CXX) $(CXXFLAGS) $< -o $@ $(LFLAGS)

mandelbrot: mandelbrot.cpp
	$(CXX) $(CXXFLAGS) $< -o $@ $(LFLAGS)

nqueens: nqueens.cpp
	$(CXX) $(CXXFLAGS) $< -o $@ $(LFLAGS)

sierpinski: sierpinski.cpp
	$(CXX) $(CXXFLAGS) $< -o $@ $(LFLAGS)

tree: tree.cpp
	$(CXX) $(CXXFLAGS) $< -o $@ $(LFLAGS)

cantor: cantor.cpp
	$(CXX) $(CXXFLAGS) $< -o $@ $(LFLAGS)

hilbert: hilbert.cpp
	$(CXX) $(CXXFLAGS) $< -o $@ $(LFLAGS)
