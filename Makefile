CXX=g++
CXXFLAGS=-std=c++11 -pedantic -Wall -Wextra -Werror -static
CXXOPTFLAGS=-O2 -fno-rtti

arps.o: arps.cpp arps.hpp
	$(CXX) $(CXXFLAGS) $(CXXOPTFLAGS) -c arps.cpp

bestfit.o: bestfit.cpp bestfit.hpp
	$(CXX) $(CXXFLAGS) $(CXXOPTFLAGS) -c bestfit.cpp

convex.o: convex.cpp convex.hpp
	$(CXX) $(CXXFLAGS) $(CXXOPTFLAGS) -c convex.cpp
