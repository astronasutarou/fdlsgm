CC     := g++
TESTPY := --repository-url https://test.pypi.org/legacy/
LIBS   := -I./src -L./src -lm
CFLAGS := -std=c++11 -g -Wall -O3
CXX    := $(CC) $(LIBS) $(CFLAGS)

HEADER :=
SOURCE := $(wildcard src/*.cc)
OBJECT := $(patsubst %.cc,%.o,$(SOURCE))

.PHONY: clean build build_pypi upload_test upload_pypi

all:

.cc.o: $(HEADER)
	$(CXX) -o $@ -c $<

test/%:test/%.cc $(OBJECT)
	$(CXX) -o $@ $^

build:
	python setup.py build_ext --inplace

build_pypi: build
	python setup.py sdist bdist_wheel -p manylinux1_x86_64

upload_test: build_pypi
	twine upload --skip-existing $(TESTPY) dist/*

upload_pypi: build_pypi
	twine upload --skip-existing dist/*

clean:
	rm -r $(OBJECT)
