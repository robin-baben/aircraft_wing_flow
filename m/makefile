run: main
	./main
main: main.cpp point3D.o frame.o
	g++ main.cpp point3D.o frame.o -lblas -llapack -fopenmp -O3 -o main
point3D.o: point3D.cpp point3D.h
	g++ -c point3D.cpp -o point3D.o
frame.o: frame.cpp frame.h
	g++ -c frame.cpp -o frame.o