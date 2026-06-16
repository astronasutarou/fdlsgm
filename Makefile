CC     := g++
TESTPY := --repository-url https://test.pypi.org/legacy/
LIBS   := -I./src -L./src -lm
CFLAGS := -std=c++11 -g -Wall -O3
CXX    := $(CC) $(LIBS) $(CFLAGS)

HEADER :=
SOURCE := $(wildcard src/*.cc)
OBJECT := $(patsubst %.cc,%.o,$(SOURCE))
EGG_INFO := fdlsgm.egg-info
BUILD_DIR := build
BUILD_LIB := $(BUILD_DIR)/lib
BUILD_TEMP := $(BUILD_DIR)/temp
DIST_DIR := $(BUILD_DIR)/dist

.PHONY: clean build build_pypi upload_test upload_pypi

all:

test/%: test/%.cc $(OBJECT)
	$(CXX) -o $@ $^

%.o: %.cc $(HEADER)
	$(CXX) -o $@ -c $<

build:
	python setup.py build_ext \
		--build-lib $(BUILD_LIB) --build-temp $(BUILD_TEMP)
	python -m build --outdir $(DIST_DIR)

clean:
	rm -rf $(BUILD_DIR) $(EGG_INFO) $(OBJECT)
