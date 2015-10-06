PREFIX := /usr/local

CXXSTD := c++11

ifndef CXX
  CXX := g++
endif

#CPPFLAGS  :=

ifndef CXXFLAGS
  CXXFLAGS := -O3
  ifdef DEBUG
	ifneq ($(DEBUG),0)
	  CXXFLAGS += -g
	endif
  endif
endif

ifndef PYTHON
  PYTHON := python
endif

ifndef CYTHON
  CYTHON := cython
endif

HAVE_ROOT := $(shell which root-config 2> /dev/null)

LIBHEADERS := $(wildcard include/Professor/*.h)
LIBSOURCES := $(wildcard src/*.cc)
LIBOBJECTS := $(patsubst %,obj/%.o, ParamPoints Ipol)
TESTPROGS  := test/testParamPoints test/testIpol
CYTHONSOURCES := $(wildcard pyext/professor2/*.pxd) $(wildcard pyext/professor2/*.pyx)

.PHONY := all lib pyext tests cxxtests pytests check icheck clean root


all: lib pyext tests
	@true

lib: lib/libProfessor2.so
	@true

lib/libProfessor2.so: $(LIBOBJECTS)
	@true
	$(CXX) -shared -Wl,-soname,libProfessor2.so -o $@ $(LIBOBJECTS)

obj/%.o: src/%.cc $(LIBHEADERS)
	mkdir -p obj lib
	$(CXX) -std=$(CXXSTD) -Iinclude $(CPPFLAGS) $(CXXFLAGS) -c -fPIC $< -o $@

pyext: pyext/professor2/core.so $(wildcard pyext/professor2/*.py)
	$(PYTHON) pyext/setup.py install --prefix=.

pyext/professor2/core.so: $(LIBHEADERS) $(CYTHONSOURCES) lib
	$(CYTHON) pyext/professor2/core.pyx --cplus
	$(PYTHON) pyext/setup.py build_ext -i --force

tests: cxxtests pytests
	@true

cxxtests: $(TESTPROGS)
	@true

test/%: test/%.cc $(LIBHEADERS) lib
	$(CXX) -std=$(CXXSTD) -Iinclude $(CPPFLAGS) $(CXXFLAGS) $< -Llib -lProfessor2 -o $@

ifdef $(HAVE_ROOT)
root: src/testRoot.cc $(LIBHEADERS) lib
	$(CXX) -std=$(CXXSTD) $(CPPFLAGS) $(CXXFLAGS) $< -Iinclude `root-config --cflags --libs` -Llib -lProfessor2 -o test/test$@
endif

pytests: pyext
	@true

check: tests
	@echo
	@echo "testParamPoints" && test/testParamPoints && echo "\n\n"
	@echo "testIpol" && test/testIpol && echo "\n\n"

icheck: tests
	test/testPython
	test/testPython1D
	test/testPython2D

clean:
	rm -rf obj/*.o lib/*
	rm -f pyext/professor2/core.cpp pyext/professor2/core.so
	rm -f $(TESTPROGS)

install: all
	mkdir -p $(PREFIX)/bin && cp bin/* $(PREFIX)/bin/
	mkdir -p $(PREFIX)/include && cp -r include/Professor $(PREFIX)/include/
	test -d lib   && mkdir -p $(PREFIX)/lib   && cp -r lib/* $(PREFIX)/lib/ || true
	test -d lib64 && mkdir -p $(PREFIX)/lib64 && cp -r lib64/* $(PREFIX)/lib64/ || true
	cp setup.sh $(PREFIX)
