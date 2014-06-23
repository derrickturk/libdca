CXX=g++
CXXFLAGS=-std=c++11 -pedantic -Wall -Wextra -Werror -static
CXXOPTFLAGS=-O2 -fno-rtti -msse3 -mfpmath=sse -ffast-math
CONFIG=-DDCA_IOSTREAMS

typecurve.exe: typecurve.cpp *.hpp
	$(CXX) $(CXXFLAGS) $(CONFIG) $(CXXOPTFLAGS) -o typecurve typecurve.cpp

peakmonth.exe: peakmonth.cpp *.hpp
	$(CXX) $(CXXFLAGS) $(CONFIG) $(CXXOPTFLAGS) -o peakmonth peakmonth.cpp

fit.exe: fit.cpp decline.hpp exponential.hpp hyperbolic.hpp hyptoexp.hpp bestfit.hpp production.hpp
	$(CXX) $(CXXFLAGS) $(CONFIG) $(CXXOPTFLAGS) -o fit fit.cpp

arps.exe: arps.cpp decline.hpp exponential.hpp hyperbolic.hpp hyptoexp.hpp bestfit.hpp any_decline.hpp
	$(CXX) $(CXXFLAGS) $(CONFIG) $(CXXOPTFLAGS) -o arps arps.cpp

bestfit.o: bestfit.cpp bestfit.hpp
	$(CXX) $(CXXFLAGS) $(CONFIG) $(CXXOPTFLAGS) -c bestfit.cpp

convex.exe: convex.cpp convex.hpp
	$(CXX) $(CXXFLAGS) $(CONFIG) $(CXXOPTFLAGS) -o convex convex.cpp

production.exe: production.cpp production.hpp
	$(CXX) $(CXXFLAGS) $(CONFIG) $(CXXOPTFLAGS) -o production production.cpp

clean:
	-rm *.exe *.o *.a
