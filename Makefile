all: a.out

OBJECTS = task1.o reader.o Gauss.o matrix.o lu.o

%.o : %.cpp
	g++ -O3 -std=c++11 -c $<

a.out: $(OBJECTS)
	g++ -o a.out $(OBJECTS)