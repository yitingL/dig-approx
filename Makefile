CXX = g++
CXXFLAGS = -Wall -std=c++17 -w -O3 -lz

COMMON_SRC = $(wildcard src/common/*.cpp)
COMMON_OBJ = $(COMMON_SRC:.cpp=.o)
ALS_SRC = $(wildcard src/als/*.cpp)
ALS_OBJ = $(ALS_SRC:.cpp=.o)
COMP_SRC = $(wildcard src/comp/*.cpp)
COMP_OBJ = $(COMP_SRC:.cpp=.o)

TARGET = main

all: main
	chmod 777 ext/abc
	chmod 777 ext/tig
main : $(ALS_OBJ) $(COMMON_OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<
clean:
	rm -f main comp $(ALS_OBJ) $(COMP_OBJ) $(COMMON_OBJ)