.PHONY: all clean strip

cfiles = optimize.h utils.h config.h Makefile

rkhfiles = rk.h \
	rk_2_2.h \
	rk_3_2.h rk_3_3.h \
	rk_4_2.h rk_4_3.h rk_4_4.h \
	rk_5_2.h rk_5_3.h rk_5_4.h \
	rk_6_2.h rk_6_3.h rk_6_4.h

hfiles = steps.h $(rkhfiles)

rkcfiles = rk.c \
	rk_2_2.c \
	rk_3_2.c rk_3_3.c \
	rk_4_2.c rk_4_3.c rk_4_4.c \
	rk_5_2.c rk_5_3.c rk_5_4.c \
	rk_6_2.c rk_6_3.c rk_6_4.c

rkofiles = rk.o \
	rk_2_2.o \
	rk_3_2.o rk_3_3.o \
	rk_4_2.o rk_4_3.o rk_4_4.o \
	rk_5_2.o rk_5_3.o rk_5_4.o \
	rk_6_2.o rk_6_3.o rk_6_4.o

ofiles = optimize.o steps.o $(rkofiles) utils.o

cc = gcc -flto -g -std=gnu11
cflags = `pkg-config --cflags libxml-2.0 glib-2.0 gthread-2.0 gsl` -c -O3 \
	-march=native -Wall -Wextra -D_FORTIFY_SOURCE=2
libs = -lm `pkg-config --libs libxml-2.0 glib-2.0 gthread-2.0 gsl`

all: ode ode.pdf

ode: $(ofiles) ode.o
	$(cc) ode.o $(ofiles) $(libs) -o ode

utils.o: utils.c utils.h config.h
	$(cc) $(cflags) utils.c -o utils.o

optimize.o: optimize.c optimize.h $(cfiles)
	$(cc) $(cflags) optimize.c -o optimize.o

steps.o: steps.c $(shfiles) $(cfiles)
	$(cc) $(cflags) steps.c -o steps.o

rk.o: rk.c $(rkhfiles) $(cfiles)
	$(cc) $(cflags) rk.c -o rk.o

rk_2_2.o: rk_2_2.c rk_2_2.h rk.h $(cfiles)
	$(cc) $(cflags) rk_2_2.c -o rk_2_2.o

rk_3_2.o: rk_3_2.c rk_3_2.h rk.h $(cfiles)
	$(cc) $(cflags) rk_3_2.c -o rk_3_2.o

rk_3_3.o: rk_3_3.c rk_3_3.h rk_3_2.h rk.h $(cfiles)
	$(cc) $(cflags) rk_3_3.c -o rk_3_3.o

rk_4_2.o: rk_4_2.c rk_4_2.h rk.h $(cfiles)
	$(cc) $(cflags) rk_4_2.c -o rk_4_2.o

rk_4_3.o: rk_4_3.c rk_4_3.h rk_4_2.h rk.h $(cfiles)
	$(cc) $(cflags) rk_4_3.c -o rk_4_3.o

rk_4_4.o: rk_4_4.c rk_4_4.h rk_4_3.h rk.h $(cfiles)
	$(cc) $(cflags) rk_4_4.c -o rk_4_4.o

rk_5_2.o: rk_5_2.c rk_5_2.h rk.h $(cfiles)
	$(cc) $(cflags) rk_5_2.c -o rk_5_2.o

rk_5_3.o: rk_5_3.c rk_5_3.h rk.h $(cfiles)
	$(cc) $(cflags) rk_5_3.c -o rk_5_3.o

rk_5_4.o: rk_5_4.c rk_5_4.h rk.h $(cfiles)
	$(cc) $(cflags) rk_5_4.c -o rk_5_4.o

rk_6_2.o: rk_6_2.c rk_6_2.h rk.h $(cfiles)
	$(cc) $(cflags) rk_6_2.c -o rk_6_2.o

rk_6_3.o: rk_6_3.c rk_6_3.h rk.h $(cfiles)
	$(cc) $(cflags) rk_6_3.c -o rk_6_3.o

rk_6_4.o: rk_6_4.c rk_6_4.h rk.h $(cfiles)
	$(cc) $(cflags) rk_6_4.c -o rk_6_4.o

ode.o: ode.c $(hfiles) $(cfiles)
	$(cc) $(cflags) ode.c -o ode.o

ode.pdf: ode.tex Makefile
	pdflatex ode
	pdflatex ode
	pdflatex ode

clean:
	rm -rf ode *.{o,aux,toc,log} html latex

strip: ode
	strip ode
