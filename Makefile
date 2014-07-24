CXX?=g++

CXXOPTFLAGS=-O2 -msse3 -mfpmath=sse -ffast-math
RTTI?=-fno-rtti
CXXFLAGS=-std=c++11 -pedantic -Wall -Wextra -Werror $(CXXOPTFLAGS)

INCLUDEDIR=include
LDFLAGS=-static

TESTINCLUDE=
TESTLIB=-Wl,-Bstatic -lboost_unit_test_framework

CONFIG=
#CONFIG=-DDCA_NO_IOSTREAMS

# Mac OS X is a nightmarish hellscape, so hold on to your butts

ifeq ($(OS),Windows_NT)
    # we're hosed anyway
else
    UNAME_S := $(shell uname -s)
    ifeq ($(UNAME_S),Darwin)
        LDFLAGS=
        TESTINCLUDE=/opt/local/include
        TESTLIB=/opt/local/lib/libboost_unit_test_framework-mt.a
    endif
endif

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

TESTS := $(patsubst %.cpp,%,$(wildcard tests/*.cpp))

examples: $(EXAMPLES)

$(EXAMPLES): %: %.cpp $(INCLUDES)
	$(CXX) -I$(INCLUDEDIR) $(CXXFLAGS) $(RTTI) -o $@ $< $(LDFLAGS)

tests: $(TESTS)

# rtti mandatory for boost.test
$(TESTS): %: %.cpp $(INCLUDES)
	$(CXX) -I$(INCLUDEDIR) -I$(TESTINCLUDE) $(CXXFLAGS) -o $@ $< $(LDFLAGS) $(TESTLIB)

clean:
	-rm *.o *.a
