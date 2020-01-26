CC     := g++
TWINE  := twine upload --skip-existing
LIBS   := -I./src -L./src -lm
CFLAGS := -std=c++11 -g -Wall -O3
CXX    := $(CC) $(LIBS) $(CFLAGS)

HEADER := src/fdlsgm.h src/linalg.h
SOURCE := $(wildcard src/*.cc)
OBJECT := $(patsubst %.cc,%.o,$(SOURCE))

.PHONY: clean build_pypi upload_pypi

all: test/sample_xyt test/sample_xyz

test/sample_%: test/sample_%.cc $(OBJECT) $(HEADER)
	$(CXX) -o $@ $< $(OBJECT)

.cc.o: $(HEADER)
	$(CXX) -o $@ -c $<

build_pypi:
	python setup.py sdist bdist_wheel -p manylinux1_x86_64

upload_test:
	$(TWINE) --repository-url https://test.pypi.org/legacy/ dist/*

upload_pypi:
	$(TWINE) --skip-existing dist/*

clean:
	rm -r $(OBJECT)
