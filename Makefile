CC     := g++
LIBS   := -I./src -L./src -lm
CFLAGS := -std=c++11 -g -Wall -O3
CXX    := $(CC) $(LIBS) $(CFLAGS)

HEADER := src/fdlsgm.h src/linalg.h
SOURCE := $(wildcard src/*.cc)
OBJECT := $(patsubst %.cc,%.o,$(SOURCE))

.PHONY: clean

all: test/sample

test/sample: test/sample.cc $(OBJECT) $(HEADER)
	$(CXX) -o $@ $< $(OBJECT)

.cc.o: $(HEADER)
	$(CXX) -o $@ -c $<

clean:
	rm -r $(OBJECT)
