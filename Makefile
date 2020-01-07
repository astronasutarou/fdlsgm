CC     := g++
LIBS   := -I./src -L./src -lm
CFLAGS := -std=c++11 -g -Wall
CXX    := $(CC) $(LIBS) $(CFLAGS)

HEADER := src/fdlsgm.h
SOURCE := $(wildcard src/*.cc)
OBJECT := $(patsubst %.cc,%.o,$(SOURCE))

.PHONY: clean test

all: test

test: test.cc $(OBJECT) $(HEADER)
	$(CXX) -o $@ $< $(OBJECT)

.cc.o: $(HEADER)
	$(CXX) -o $@ -c $<

clean:
	rm -r $(OBJECT)
