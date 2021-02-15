# MAkefile

CC = g++
CFLAGS = -O3 -Wall
ITPPFLAG = -litpp

.PHONY: all
all: mimo siso

siso: main_siso.o mymodulator.o
	$(CC) $(CFLAGS) -o $@ main_siso.o mymodulator.o $(ITPPFLAG)

mimo: main_mimo.o mymodulator.o
	$(CC) $(CFLAGS) -o $@ main_mimo.o mymodulator.o $(ITPPFLAG)

main_mimo.o: main_mimo.cpp
	$(CC) $(CFLAGS) -c main_mimo.cpp $(ITPPFLAG)

main_siso.o: main_siso.cpp
	$(CC) $(CFLAGS) -c main_siso.cpp $(ITPPFLAG)

mymodulator.o: mymodulator.cpp mymodulator.h
	$(CC) $(CFLAGS) -o $@ -c mymodulator.cpp $(ITPPFLAG)

.PHONY: clean
clean:
	rm -f mimo siso *.o
