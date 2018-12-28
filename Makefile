.PHONY: all clean strip

cfiles = optimize.h utils.h config.h Makefile

shfiles = steps.h \
	steps_3_2.h steps_3_3.h \
	steps_4_2.h steps_4_3.h \
	steps_5_2.h steps_5_3.h steps_5_4.h \
	steps_6_2.h steps_6_3.h steps_6_4.h steps_6_5.h \
	steps_7_2.h steps_7_3.h steps_7_4.h steps_7_5.h \
	steps_8_2.h steps_8_3.h steps_8_4.h steps_8_5.h steps_8_6.h \
	steps_9_2.h steps_9_3.h steps_9_4.h steps_9_5.h steps_9_6.h \
	steps_10_2.h steps_10_3.h steps_10_4.h steps_10_5.h steps_10_6.h \
	steps_11_2.h steps_11_3.h steps_11_4.h steps_11_5.h steps_11_6.h \

rkhfiles = rk.h \
	rk_2_2.h \
	rk_3_2.h rk_3_3.h \
	rk_4_2.h rk_4_3.h rk_4_4.h \
	rk_5_2.h rk_5_3.h rk_5_4.h \
	rk_6_2.h rk_6_3.h rk_6_4.h

hfiles = $(shfiles) $(rkhfiles)

scfiles = steps.c \
	steps_3_2.c steps_3_3.c \
	steps_4_2.c steps_4_3.c \
	steps_5_2.c steps_5_3.c steps_5_4.c \
	steps_6_2.c steps_6_3.c steps_6_4.c steps_6_5.c \
	steps_7_2.c steps_7_3.c steps_7_4.c steps_7_5.c \
	steps_8_2.c steps_8_3.c steps_8_4.c steps_8_5.c steps_8_6.c \
	steps_9_2.c steps_9_3.c steps_9_4.c steps_9_5.c steps_9_6.c \
	steps_10_2.c steps_10_3.c steps_10_4.c steps_10_5.c steps_10_6.c \
	steps_11_2.c steps_11_3.c steps_11_4.c steps_11_5.c steps_11_6.c \

rkcfiles = rk.c \
	rk_2_2.c \
	rk_3_2.c rk_3_3.c \
	rk_4_2.c rk_4_3.c rk_4_4.c \
	rk_5_2.c rk_5_3.c rk_5_4.c \
	rk_6_2.c rk_6_3.c rk_6_4.c

sofiles = steps.o \
	steps_3_2.o steps_3_3.o \
	steps_4_2.o steps_4_3.o \
	steps_5_2.o steps_5_3.o steps_5_4.o \
	steps_6_2.o steps_6_3.o steps_6_4.o steps_6_5.o \
	steps_7_2.o steps_7_3.o steps_7_4.o steps_7_5.o \
	steps_8_2.o steps_8_3.o steps_8_4.o steps_8_5.o steps_8_6.o \
	steps_9_2.o steps_9_3.o steps_9_4.o steps_9_5.o steps_9_6.o \
	steps_10_2.o steps_10_3.o steps_10_4.o steps_10_5.o steps_10_6.o \
	steps_11_2.o steps_11_3.o steps_11_4.o steps_11_5.o steps_11_6.o \

rkofiles = rk.o \
	rk_2_2.o \
	rk_3_2.o rk_3_3.o \
	rk_4_2.o rk_4_3.o rk_4_4.o \
	rk_5_2.o rk_5_3.o rk_5_4.o \
	rk_6_2.o rk_6_3.o rk_6_4.o

ofiles = optimize.o $(sofiles) $(rkofiles) utils.o

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

steps_3_2.o: steps_3_2.c steps_3_2.h steps.h $(cfiles)
	$(cc) $(cflags) steps_3_2.c -o steps_3_2.o

steps_3_3.o: steps_3_3.c steps_3_3.h steps_3_2.h steps.h $(cfiles)
	$(cc) $(cflags) steps_3_3.c -o steps_3_3.o

steps_4_2.o: steps_4_2.c steps_4_2.h steps.h $(cfiles)
	$(cc) $(cflags) steps_4_2.c -o steps_4_2.o

steps_4_3.o: steps_4_3.c steps_4_3.h steps_4_2.h steps.h $(cfiles)
	$(cc) $(cflags) steps_4_3.c -o steps_4_3.o

steps_5_2.o: steps_5_2.c steps_5_2.h steps.h $(cfiles)
	$(cc) $(cflags) steps_5_2.c -o steps_5_2.o

steps_5_3.o: steps_5_3.c steps_5_3.h steps_5_2.h steps.h $(cfiles)
	$(cc) $(cflags) steps_5_3.c -o steps_5_3.o

steps_5_4.o: steps_5_4.c steps_5_4.h steps_5_3.h steps.h $(cfiles)
	$(cc) $(cflags) steps_5_4.c -o steps_5_4.o

steps_6_2.o: steps_6_2.c steps_6_2.h steps.h $(cfiles)
	$(cc) $(cflags) steps_6_2.c -o steps_6_2.o

steps_6_3.o: steps_6_3.c steps_6_3.h steps_6_2.h steps.h $(cfiles)
	$(cc) $(cflags) steps_6_3.c -o steps_6_3.o

steps_6_4.o: steps_6_4.c steps_6_4.h steps_6_3.h steps.h $(cfiles)
	$(cc) $(cflags) steps_6_4.c -o steps_6_4.o

steps_6_5.o: steps_6_5.c steps_6_5.h steps_6_4.h steps.h $(cfiles)
	$(cc) $(cflags) steps_6_5.c -o steps_6_5.o

steps_7_2.o: steps_7_2.c steps_7_2.h steps.h $(cfiles)
	$(cc) $(cflags) steps_7_2.c -o steps_7_2.o

steps_7_3.o: steps_7_3.c steps_7_3.h steps_7_2.h steps.h $(cfiles)
	$(cc) $(cflags) steps_7_3.c -o steps_7_3.o

steps_7_4.o: steps_7_4.c steps_7_4.h steps_7_3.h steps.h $(cfiles)
	$(cc) $(cflags) steps_7_4.c -o steps_7_4.o

steps_7_5.o: steps_7_5.c steps_7_5.h steps_7_4.h steps.h $(cfiles)
	$(cc) $(cflags) steps_7_5.c -o steps_7_5.o

steps_8_2.o: steps_8_2.c steps_8_2.h steps.h $(cfiles)
	$(cc) $(cflags) steps_8_2.c -o steps_8_2.o

steps_8_3.o: steps_8_3.c steps_8_3.h steps_8_2.h steps.h $(cfiles)
	$(cc) $(cflags) steps_8_3.c -o steps_8_3.o

steps_8_4.o: steps_8_4.c steps_8_4.h steps_8_3.h steps.h $(cfiles)
	$(cc) $(cflags) steps_8_4.c -o steps_8_4.o

steps_8_5.o: steps_8_5.c steps_8_5.h steps_8_4.h steps.h $(cfiles)
	$(cc) $(cflags) steps_8_5.c -o steps_8_5.o

steps_8_6.o: steps_8_6.c steps_8_6.h steps_8_5.h steps.h $(cfiles)
	$(cc) $(cflags) steps_8_6.c -o steps_8_6.o

steps_9_2.o: steps_9_2.c steps_9_2.h steps.h $(cfiles)
	$(cc) $(cflags) steps_9_2.c -o steps_9_2.o

steps_9_3.o: steps_9_3.c steps_9_3.h steps_9_2.h steps.h $(cfiles)
	$(cc) $(cflags) steps_9_3.c -o steps_9_3.o

steps_9_4.o: steps_9_4.c steps_9_4.h steps_9_3.h steps.h $(cfiles)
	$(cc) $(cflags) steps_9_4.c -o steps_9_4.o

steps_9_5.o: steps_9_5.c steps_9_5.h steps_9_4.h steps.h $(cfiles)
	$(cc) $(cflags) steps_9_5.c -o steps_9_5.o

steps_9_6.o: steps_9_6.c steps_9_6.h steps_9_5.h steps.h $(cfiles)
	$(cc) $(cflags) steps_9_6.c -o steps_9_6.o

steps_10_2.o: steps_10_2.c steps_10_2.h steps.h $(cfiles)
	$(cc) $(cflags) steps_10_2.c -o steps_10_2.o

steps_10_3.o: steps_10_3.c steps_10_3.h steps_10_2.h steps.h $(cfiles)
	$(cc) $(cflags) steps_10_3.c -o steps_10_3.o

steps_10_4.o: steps_10_4.c steps_10_4.h steps_10_3.h steps.h $(cfiles)
	$(cc) $(cflags) steps_10_4.c -o steps_10_4.o

steps_10_5.o: steps_10_5.c steps_10_5.h steps_10_4.h steps.h $(cfiles)
	$(cc) $(cflags) steps_10_5.c -o steps_10_5.o

steps_10_6.o: steps_10_6.c steps_10_6.h steps_10_5.h steps.h $(cfiles)
	$(cc) $(cflags) steps_10_6.c -o steps_10_6.o

steps_11_2.o: steps_11_2.c steps_11_2.h steps.h $(cfiles)
	$(cc) $(cflags) steps_11_2.c -o steps_11_2.o

steps_11_3.o: steps_11_3.c steps_11_3.h steps_11_2.h steps.h $(cfiles)
	$(cc) $(cflags) steps_11_3.c -o steps_11_3.o

steps_11_4.o: steps_11_4.c steps_11_4.h steps_11_3.h steps.h $(cfiles)
	$(cc) $(cflags) steps_11_4.c -o steps_11_4.o

steps_11_5.o: steps_11_5.c steps_11_5.h steps_11_4.h steps.h $(cfiles)
	$(cc) $(cflags) steps_11_5.c -o steps_11_5.o

steps_11_6.o: steps_11_6.c steps_11_6.h steps_11_5.h steps.h $(cfiles)
	$(cc) $(cflags) steps_11_6.c -o steps_11_6.o

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
