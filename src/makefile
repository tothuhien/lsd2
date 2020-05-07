CXX=g++
CXXFLAGS=-std=c++11 -Ofast
EXEC=lsd2
SRC= $(wildcard *.cpp)
OBJ= $(SRC:.cpp=.o)

all: $(EXEC)

$(EXEC): $(OBJ)
	$(CXX) -o $@ $^ 

%.o: %.cpp %.h
	$(CXX) $(CXXFLAGS) -o $@ -c $< 
	
.PHONY: clean mrproper

clean:
	rm -rf *.o

mrproper: clean
	rm -rf $(EXEC)
