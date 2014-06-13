CXX=g++
CXXFLAGS=-std=c++11 -pedantic -Wall -Wextra -Werror -static
CXXOPTFLAGS=-O2 -fno-rtti -msse3 -mfpmath=sse -ffast-math

arps.exe: arps.cpp decline.hpp exponential.hpp
	$(CXX) $(CXXFLAGS) $(CXXOPTFLAGS) -o arps arps.cpp

bestfit.o: bestfit.cpp bestfit.hpp
	$(CXX) $(CXXFLAGS) $(CXXOPTFLAGS) -c bestfit.cpp

convex.exe: convex.cpp convex.hpp
	$(CXX) $(CXXFLAGS) $(CXXOPTFLAGS) -o convex convex.cpp
