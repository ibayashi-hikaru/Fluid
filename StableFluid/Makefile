SRCS=$(shell find ./ -name "*.cpp")
OBJS=$(SRCS:.cpp=.o) #Foo.o
DEPS=$(SRCS:.cpp=.d) #Foo.d

CXX=clang++
CXXFLAGS=-O0 -g -std=c++11 -Wno-deprecated -framework GLUT -framework OpenGL -I. -I/usr/local/include -MMD # -O4
ASMFLAGS=-c
LDFLAGS=-lm -lobjc -L. -L/usr/X11R6/lib -L/usr/lib
TARGET=Fluid

.PHONY: all clean
all: $(TARGET)

$(TARGET): $(OBJS) # Classic: Foo.o
	$(CXX) $(OBJS) -o $(TARGET) $(CXXFLAGS) $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(ASMFLAGS) $(CXXFLAGS) $<

clean:
	rm -rf $(OBJS) $(TARGET) $(DEPS)

include $(shell find ./ -name "*.d") # Foo.o: Foo.cpp Bar.o
