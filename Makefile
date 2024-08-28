# Makefile for disccofan project

# Compilers
CC = gcc
CXX = g++
MPICC = mpicc
MPICXX = mpicxx

# Compiler Flags
CFLAGS = -O3 -fopenmp -pedantic -Wall
CXXFLAGS = -std=c++11 -O3 -fopenmp -pedantic -Wall

# Directories
SRC_DIR = src/sources
INC_DIR = src/headers
OBJ_DIR = obj

# Libraries
LIBS = -fopenmp -lfreeimage -lhdf5 -lcfitsio -lconfig -lm
MPI_LIBS = $(shell $(MPICXX) -showme:link)

# Include directories
INCLUDES = -I$(INC_DIR) -I/usr/include/cfitsio -I/usr/include
MPI_INCLUDES = $(shell $(MPICXX) -showme:compile)

# Source and Object files
SOURCES = $(wildcard $(SRC_DIR)/*.c)
OBJECTS = $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(SOURCES))

# Executable
TARGET = disccofan

# Create directories if they do not exist
$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)


# Build target
$(TARGET): $(OBJ_DIR) $(OBJECTS)
	$(MPICC) $(OBJECTS) -o $(TARGET) $(LIBS) $(MPI_LIBS)

# Build object files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c | $(OBJ_DIR)
	$(MPICC) $(CFLAGS) $(INCLUDES) $(MPI_INCLUDES) -c $< -o $@

# Clean
clean:
	rm -rf $(OBJ_DIR) $(TARGET)

.PHONY: clean

# Debugging information
print-%: ; @echo $* = $($*)
