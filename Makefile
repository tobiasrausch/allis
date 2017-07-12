DEBUG ?= 0
STATIC ?= 0

# Submodules
PWD = $(shell pwd)
EBROOTHTSLIB ?= ${PWD}/src/htslib/
BOOST_ROOT ?= ${PWD}/src/modular-boost/

# Install dir
prefix = ${PWD}
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin

# Flags
CXX=g++
CXXFLAGS += -isystem ${EBROOTHTSLIB} -isystem ${BOOST_ROOT} -pedantic -W -Wall -Wno-unknown-pragmas -D__STDC_LIMIT_MACROS -fno-strict-aliasing
LDFLAGS += -L${EBROOTHTSLIB} -L${EBROOTHTSLIB}/lib -L${BOOST_ROOT}/stage/lib -lboost_iostreams -lboost_filesystem -lboost_system -lboost_program_options -lboost_date_time

# Additional flags for release/debug
ifeq (${STATIC}, 1)
	LDFLAGS += -static -static-libgcc -pthread -lhts -lz -llzma -lbz2
else
	LDFLAGS += -lhts -lz -llzma -lbz2 -Wl,-rpath,${EBROOTHTSLIB},-rpath,${EBROOTHTSLIB}/lib,-rpath,${BOOST_ROOT}/stage/lib
endif
ifeq (${DEBUG}, 1)
	CXXFLAGS += -g -O0 -fno-inline -DDEBUG
else ifeq (${DEBUG}, 2)
	CXXFLAGS += -g -O0 -fno-inline -DPROFILE
	LDFLAGS += -lprofiler -ltcmalloc
else
	CXXFLAGS += -O3 -fno-tree-vectorize -DNDEBUG
endif
ifeq (${EBROOTHTSLIB}, ${PWD}/src/htslib/)
	SUBMODULES += .htslib
endif
ifeq (${BOOST_ROOT}, ${PWD}/src/modular-boost/)
	SUBMODULES += .boost
endif


# External sources
BOOSTSOURCES = $(wildcard src/modular-boost/libs/iostreams/include/boost/iostreams/*.hpp)
HTSLIBSOURCES = $(wildcard src/htslib/*.c) $(wildcard src/htslib/*.h)
SOURCES = $(wildcard src/*.h) $(wildcard src/*.cpp)

# Targets
BUILT_PROGRAMS = src/allis
TARGETS = ${SUBMODULES} ${BUILT_PROGRAMS}

all:   	$(TARGETS)

.htslib: $(HTSLIBSOURCES)
	cd src/htslib && make && make lib-static && cd ../../ && touch .htslib

.boost: $(BOOSTSOURCES)
	cd src/modular-boost && ./bootstrap.sh --prefix=${PWD}/src/modular-boost --without-icu --with-libraries=iostreams,filesystem,system,program_options,date_time && ./b2 && ./b2 headers && cd ../../ && touch .boost

src/allis: ${SUBMODULES} $(SOURCES)
	$(CXX) $(CXXFLAGS) $@.cpp -o $@ $(LDFLAGS)

install: ${BUILT_PROGRAMS}
	mkdir -p ${bindir}
	install -p ${BUILT_PROGRAMS} ${bindir}

clean:
	cd src/htslib && make clean
	cd src/modular-boost && ./b2 --clean-all
	rm -f $(TARGETS) $(TARGETS:=.o)