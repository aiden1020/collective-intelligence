CXX = g++
CXXFLAGS = -Wall -Wextra -std=c++11

TARGET = PSO

SRCS = main.cpp

all: $(TARGET)

$(TARGET): $(SRCS)
	$(CXX) $(CXXFLAGS) -o $@ $^

clean:
	rm -f $(TARGET)
