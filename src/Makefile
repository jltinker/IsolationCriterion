LIB = -lm


CC = gcc
CFLAGS = -O2

OBJS10 = isolation_condition.o file_control.o nrutil.o midpnt.o qromo.o polint.o
isolation_condition:	$(OBJS10)
	$(CC) -o $@ $(OBJS10) $(LIB)
#	cp -f $@ $(HOME)/exec/$@

clean:
	rm -f *.o
