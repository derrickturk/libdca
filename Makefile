CXX?=g++

CXXOPTFLAGS=-O2 -msse3 -mfpmath=sse -ffast-math
RTTI?=-fno-rtti
CXXFLAGS=-std=c++14 -pedantic -Wall -Wextra -Werror $(CXXOPTFLAGS)

INCLUDEDIR=include
LDFLAGS=-static

BOOSTINCLUDE=
BOOSTTESTLINK=-Wl,-Bstatic -lboost_unit_test_framework
BOOSTFLAGS=-Wno-deprecated-declarations

CONFIG=
#CONFIG=-DDCA_NO_IOSTREAMS

# Mac OS X is a nightmarish hellscape, so hold on to your butts

ifeq ($(OS),Windows_NT)
    CLEAN=clean-win
    BOOSTINCLUDE=-IC:/Users/Derrick/local/boost/boost_1_59_0
    BOOSTTESTLINK=-LC:/Users/Derrick/local/boost/boost_1_59_0/stage/lib -Wl,-Bstatic -lboost_unit_test_framework-mgw51-mt-1_59
else
    CLEAN=clean-unix
    UNAME_S := $(shell uname -s)
    ifeq ($(UNAME_S),Darwin)
        LDFLAGS=
        BOOSTINCLUDE=-I/opt/local/include
        BOOSTTESTLINK=/opt/local/lib/libboost_unit_test_framework-mt.a
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
	$(CXX) -I$(INCLUDEDIR) $(CXXFLAGS) $(CONFIG) $(RTTI) -o $@ $< $(LDFLAGS)

tests: $(TESTS)

# rtti mandatory for boost.test
$(TESTS): %: %.cpp $(INCLUDES)
	$(CXX) -I$(INCLUDEDIR) $(BOOSTINCLUDE) $(CXXFLAGS) $(BOOSTFLAGS) $(CONFIG) -o $@ $< $(LDFLAGS) $(BOOSTTESTLINK)

clean: $(CLEAN)

clean-unix:
	-rm $(TESTS)
	-rm $(EXAMPLES)

clean-win:
	-rm tests/*.exe
	-rm examples/*.exe
