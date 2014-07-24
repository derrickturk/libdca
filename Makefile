CXX?=g++
CXXFLAGS=-std=c++11 -pedantic -Wall -Wextra -Werror
CXXOPTFLAGS=-O2 -fno-rtti -msse3 -mfpmath=sse -ffast-math
LDFLAGS=-static
INCLUDEDIR=include
CONFIG=
#CONFIG=-DDCA_NO_IOSTREAMS

.SUFFIXES:

INCLUDES=\
	$(INCLUDEDIR)/dca/any_decline.hpp \
	$(INCLUDEDIR)/dca/bestfit.hpp \
	$(INCLUDEDIR)/dca/convex.hpp \
	$(INCLUDEDIR)/dca/decline.hpp \
	$(INCLUDEDIR)/dca/exponential.hpp \
	$(INCLUDEDIR)/dca/hyperbolic.hpp \
	$(INCLUDEDIR)/dca/hyptoexp.hpp \
	$(INCLUDEDIR)/dca/production.hpp \
	$(INCLUDEDIR)/dca/tuple_tools.hpp

EXAMPLES := $(patsubst %.cpp,%,$(wildcard examples/*.cpp))

examples: $(EXAMPLES)

$(EXAMPLES): %: %.cpp $(INCLUDES)
	$(CXX) -I$(INCLUDEDIR) $(CXXFLAGS) $(CXXOPTFLAGS) -o $@ $< $(LDFLAGS)

clean:
	-rm *.o *.a
