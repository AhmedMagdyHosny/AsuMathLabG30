all: matrix

matrix: main.o matrix.o

	g++ -std=c++11 main.o matrix.o -o matrix

main.o: main.cpp
	g++ -std=c++11 -c main.cpp

matrix.o: matrix.cpp matrix.h
	g++ -std=c++11 -c matrix.cpp

clean:
	rm *.o matrix
