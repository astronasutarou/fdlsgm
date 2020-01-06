CC     := g++
LIBS   := -I./src -L./src -lm
CFLAGS := -std=c++11 -Wall
CXX    := $(CC) $(LIBS) $(CFLAGS)

HEADER := src/fdlsgm.h
SOURCE := $(wildcard src/*.cc)
TARGET := $(patsubst %.cc,%.o,$(SOURCE))

.PHONY: clean

all: $(TARGET)
	echo $(TARGET)

.cc.o: $(HEADER)
	$(CXX) -o $@ -c $<


clean:
	rm -r $(TARGET)
