# Makefile
MAIN_SRC = main.cpp

# Compiler
CXX = g++ 

# Compiler flags
CXXFLAGS = -std=c++17 -O3 -march=native -ffast-math -Wall -pthread -Iinclude 
# Debug flags
DEBUGFLAGS = -std=c++17 -g -Wall -pthread -Iinclude 

# Linker flags (libraries to link against)
LDFLAGS = -lm -ltbb 

# Target executable
TARGET = a.out

# Debug target executable
DEBUG_TARGET = debug.out

# Source directory (for non-main source files)
SRC_DIR = src

# Build directory for object files
BUILD_DIR = build

# Include directory for headers
INC_DIR = include

# Test source directory
TEST_DIR = test

# Test executable
TEST_TARGET = test.out

# Main source file (located outside src directory)
#main.cpp   # ##

# Source files (non-main source files located in src directory)
SRCS = $(wildcard $(SRC_DIR)/*.cpp)

# Test source files (located in testing directory)
TEST_SRCS = $(TEST_DIR)/test.cpp

# Object files will be placed in the build directory
MAIN_OBJ = $(BUILD_DIR)/main.o
OBJS = $(SRCS:$(SRC_DIR)/%.cpp=$(BUILD_DIR)/%.o)

# Test object files will also be placed in the build directory
TEST_OBJS = $(TEST_SRCS:$(TEST_DIR)/%.cpp=$(BUILD_DIR)/%.o)

# Default target
all: create_build_dir $(TARGET)

# Create the build directory if it doesn't exist
create_build_dir:
	@mkdir -p $(BUILD_DIR)

# Build the target executable
$(TARGET): $(MAIN_OBJ) $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(MAIN_OBJ) $(OBJS) $(LDFLAGS)

# Build the debug executable
$(DEBUG_TARGET): $(MAIN_OBJ) $(OBJS)
	$(CXX) $(DEBUGFLAGS) -o $(DEBUG_TARGET) $(MAIN_OBJ) $(OBJS) $(LDFLAGS)

# Build the test executable
$(TEST_TARGET): $(TEST_OBJS) $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TEST_TARGET) $(TEST_OBJS) $(LDFLAGS)

# Compile main.cpp into an object file
$(BUILD_DIR)/main.o: $(MAIN_SRC) $(INC_DIR)/*.h
	$(CXX) $(CXXFLAGS) -c $(MAIN_SRC) -o $(BUILD_DIR)/main.o

# Compile .cpp files from the src directory into .o files in the build directory
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp $(INC_DIR)/*.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Compile .cpp files from the testing directory into .o files in the build directory
$(BUILD_DIR)/%.o: $(TEST_DIR)/%.cpp $(INC_DIR)/*.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean build files (including the build directory)
clean:
	rm -f $(BUILD_DIR)/*.o $(TARGET) $(DEBUG_TARGET) $(TEST_TARGET)
	rmdir $(BUILD_DIR) || true

# Run the program
run: $(TARGET)
	./$(TARGET)

# Run the tests
test: create_build_dir $(TEST_TARGET)

# Debug target
debug: create_build_dir $(DEBUG_TARGET)
