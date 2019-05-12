CC=g++
CPPFLAGS=-std=c++11
EXEC=lsd2
SRC= $(wildcard *.cpp)
OBJ= $(SRC:.cpp=.o)

all: $(EXEC)

$(EXEC): $(OBJ)
	$(CC) -o $@ $^ 

%.o: %.cpp %.h
	$(CC) $(CPPFLAGS) -o $@ -c $< 
	
.PHONY: clean mrproper

clean:
	rm -rf *.o

mrproper: clean
	rm -rf $(EXEC)
