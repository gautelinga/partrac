CXX      := -h5c++
CXXFLAGS := -pedantic-errors -Wall -Wextra -Werror -O3 -std=c++17 -ffast-math -fconcepts
LDFLAGS  := -L/usr/lib -lstdc++
BUILD    := build
OBJ_DIR  := $(BUILD)
TARGET   := lbm
INCLUDE  := -Iinclude/
SRC_TRACE :=  $(wildcard src/*.cpp)

OBJECTS_TRACE := $(SRC_TRACE:%.cpp=$(OBJ_DIR)/%.o)

all: build trace

$(OBJ_DIR)/%.o: %.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ -c $<

trace: 	$(OBJECTS_TRACE)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(LDFLAGS) -o trace $(OBJECTS_TRACE)

.PHONY: all build clean allclean

build:
	@mkdir -p $(OBJ_DIR)
	ctags -e -R src/*.cpp src/*.hpp

clean:
	rm -rf build/*

allclean: clean
	rm trace
