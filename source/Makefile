CXX = g++
CXXFLAGS = -std=c++17 -Wall -O3

# Default config
CONFIG ?= K64N128

# Map CONFIG to preprocessor define
DEFS = -DCONFIG_$(CONFIG)

# Source files
SRCS = \
  argmin.cpp \
  conv.cpp \
  enc_dec.cpp \
  FeedForwardTrellis.cpp \
  lowRateListDecoder.cpp \
  minHeap.cpp \
  run_test.cpp \
  tbcc_namespace.cpp

# Object files
OBJS = $(SRCS:.cpp=.o)

# Executable name
TARGET = run_test

# Include paths (adjust if headers are in separate dirs)
INCLUDES = -I.

# Rule to build everything
all: $(TARGET)

# Compile rule
%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(DEFS) $(INCLUDES) -c $< -o $@

# Link everything
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Clean build artifacts
clean:
	rm -f $(OBJS) $(TARGET)

# Example usage: 
  # make CONFIG=K64N80
  # ./run_test -t 0
