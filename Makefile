DCFLAGS = -DDEBUG -pedantic -O3
CFLAGS = -pedantic -O3

all:
	make clean
	make fflo

fflo:	fflo.o hutils.o hgauss.o
	$(CC) $(CFLAGS)  fflo.o hutils.o hgauss.o -o fflo -lm

clean:
	rm -f *.o fflo
