CC     := g++
LIBS   := -I./src -L./src -lm
CFLAGS := -std=c++11 -g -Wall
CXX    := $(CC) $(LIBS) $(CFLAGS)

HEADER := src/fdlsgm.h
SOURCE := $(wildcard src/*.cc)
OBJECT := $(patsubst %.cc,%.o,$(SOURCE))

.PHONY: clean test

all: $(OBJECT)
	echo $(OBJECT)

test: test.cc $(OBJECT)
	$(CXX) -o $@ $^

.cc.o: $(HEADER)
	$(CXX) -o $@ -c $^


clean:
	rm -r $(OBJECT)
