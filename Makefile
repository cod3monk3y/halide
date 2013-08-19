CC=g++
CFLAGS=-I../src/include -I. `libpng-config --cflags`
LFLAGS=-L../src/bin -lHalide `libpng-config --ldflags`

test: grayscale.o
	$(CC) $(CFLAGS) $(LFLAGS) -o test grayscale.o

grayscale.o:
	$(CC) $(CFLAGS) -c grayscale.cpp

clean:
	-rm test *.o