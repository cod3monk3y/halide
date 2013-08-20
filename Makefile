CC=g++
CFLAGS=-I../src/include -I. `libpng-config --cflags`
LDFLAGS=-L../src/bin -lHalide `libpng-config --ldflags`

test: stdafx.o pockethandbook.o test.o
	$(CC) $(CFLAGS) $(LDFLAGS) -o test stdafx.o pockethandbook.o test.o

stdafx.o: stdafx.h
test.o: pockethandbook.h stdafx.h
pockethandbook.o: pockethandbook.h stdafx.h

# common header requires all .cpp files to be recompiled when the stdafx.h header changes
# %.o: stdafx.h

# for description of "$<" and "$@", see:
# www.gnu.org/software/make/manual/make.html#Automatic-Variables
# NOTE: this could be definied using predefined implicit rules
%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	-rm test *.o