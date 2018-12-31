.PHONY: all clean strip

cfiles = optimize.h utils.h config.h Makefile

rkhfiles = rk.h \
	rk_2_2.h \
	rk_3_2.h rk_3_3.h \
	rk_4_2.h rk_4_3.h rk_4_4.h \
	rk_5_2.h rk_5_3.h rk_5_4.h \
	rk_6_2.h rk_6_3.h rk_6_4.h

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

rkpgofiles = rk.pgo \
	rk_2_2.pgo \
	rk_3_2.pgo rk_3_3.pgo \
	rk_4_2.pgo rk_4_3.pgo rk_4_4.pgo \
	rk_5_2.pgo rk_5_3.pgo rk_5_4.pgo \
	rk_6_2.pgo rk_6_3.pgo rk_6_4.pgo

rkgcdafiles = rk.pgo \
	rk_2_2.pgo \
	rk_3_2.pgo rk_3_3.pgo \
	rk_4_2.pgo rk_4_3.pgo rk_4_4.pgo \
	rk_5_2.pgo rk_5_3.pgo rk_5_4.pgo \
	rk_6_2.pgo rk_6_3.pgo rk_6_4.pgo

ofiles = ode.o optimize.o steps.o $(rkofiles) utils.o

pgofiles = ode.pgo optimize.pgo steps.pgo $(rkpgofiles) utils.pgo

gcdafiles = ode.gcda optimize.gcda steps.gcda $(rkgcdafiles) utils.gcda

cc = gcc -flto -g -std=gnu11
ccgen = $(cc) -fprofile-generate
ccuse = $(cc) -fprofile-use -fprofile-correction
cflags = `pkg-config --cflags libxml-2.0 glib-2.0 gthread-2.0 gsl` -c -O3 \
	-march=native -Wall -Wextra -D_FORTIFY_SOURCE=2
libs = -lm `pkg-config --libs libxml-2.0 glib-2.0 gthread-2.0 gsl`

all: ode ode.pdf

ode: $(ofiles)
	$(ccuse) $(ofiles) $(libs) -o ode

utils.o: utils.gcda
	$(ccuse) $(cflags) utils.c -o utils.o

optimize.o: optimize.gcda
	$(ccuse) $(cflags) optimize.c -o optimize.o

steps.o: steps.gcda
	$(ccuse) $(cflags) steps.c -o steps.o

rk.o: rk.gcda
	$(ccuse) $(cflags) rk.c -o rk.o

rk_2_2.o: rk_2_2.gcda
	$(ccuse) $(cflags) rk_2_2.c -o rk_2_2.o

rk_3_2.o: rk_3_2.gcda
	$(ccuse) $(cflags) rk_3_2.c -o rk_3_2.o

rk_3_3.o: rk_3_3.gcda
	$(ccuse) $(cflags) rk_3_3.c -o rk_3_3.o

rk_4_2.o: rk_4_2.gcda
	$(ccuse) $(cflags) rk_4_2.c -o rk_4_2.o

rk_4_3.o: rk_4_3.gcda
	$(ccuse) $(cflags) rk_4_3.c -o rk_4_3.o

rk_4_4.o: rk_4_4.gcda
	$(ccuse) $(cflags) rk_4_4.c -o rk_4_4.o

rk_5_2.o: rk_5_2.gcda
	$(ccuse) $(cflags) rk_5_2.c -o rk_5_2.o

rk_5_3.o: rk_5_3.gcda
	$(ccuse) $(cflags) rk_5_3.c -o rk_5_3.o

rk_5_4.o: rk_5_4.gcda
	$(ccuse) $(cflags) rk_5_4.c -o rk_5_4.o

rk_6_2.o: rk_6_2.gcda
	$(ccuse) $(cflags) rk_6_2.c -o rk_6_2.o

rk_6_3.o: rk_6_3.gcda
	$(ccuse) $(cflags) rk_6_3.c -o rk_6_3.o

rk_6_4.o: rk_6_4.gcda
	$(ccuse) $(cflags) rk_6_4.c -o rk_6_4.o

ode.o: ode.gcda
	$(ccuse) $(cflags) ode.c -o ode.o

ode-pgo: $(pgofiles)
	$(ccgen) $(pgofiles) $(libs) -o ode-pgo

utils.pgo: utils.c utils.h config.h
	$(ccgen) $(cflags) utils.c -o utils.pgo

optimize.pgo: optimize.c optimize.h $(cfiles)
	$(ccgen) $(cflags) optimize.c -o optimize.pgo

steps.pgo: steps.c steps.h $(cfiles)
	$(ccgen) $(cflags) steps.c -o steps.pgo

rk.pgo: rk.c $(rkhfiles) $(cfiles)
	$(ccgen) $(cflags) rk.c -o rk.pgo

rk_2_2.pgo: rk_2_2.c rk_2_2.h rk.h $(cfiles)
	$(ccgen) $(cflags) rk_2_2.c -o rk_2_2.pgo

rk_3_2.pgo: rk_3_2.c rk_3_2.h rk.h $(cfiles)
	$(ccgen) $(cflags) rk_3_2.c -o rk_3_2.pgo

rk_3_3.pgo: rk_3_3.c rk_3_3.h rk_3_2.h rk.h $(cfiles)
	$(ccgen) $(cflags) rk_3_3.c -o rk_3_3.pgo

rk_4_2.pgo: rk_4_2.c rk_4_2.h rk.h $(cfiles)
	$(ccgen) $(cflags) rk_4_2.c -o rk_4_2.pgo

rk_4_3.pgo: rk_4_3.c rk_4_3.h rk_4_2.h rk.h $(cfiles)
	$(ccgen) $(cflags) rk_4_3.c -o rk_4_3.pgo

rk_4_4.pgo: rk_4_4.c rk_4_4.h rk_4_3.h rk.h $(cfiles)
	$(ccgen) $(cflags) rk_4_4.c -o rk_4_4.pgo

rk_5_2.pgo: rk_5_2.c rk_5_2.h rk.h $(cfiles)
	$(ccgen) $(cflags) rk_5_2.c -o rk_5_2.pgo

rk_5_3.pgo: rk_5_3.c rk_5_3.h rk_5_2.h rk.h $(cfiles)
	$(ccgen) $(cflags) rk_5_3.c -o rk_5_3.pgo

rk_5_4.pgo: rk_5_4.c rk_5_4.h rk_5_3.h rk.h $(cfiles)
	$(ccgen) $(cflags) rk_5_4.c -o rk_5_4.pgo

rk_6_2.pgo: rk_6_2.c rk_6_2.h rk.h $(cfiles)
	$(ccgen) $(cflags) rk_6_2.c -o rk_6_2.pgo

rk_6_3.pgo: rk_6_3.c rk_6_3.h rk_6_2.h rk.h $(cfiles)
	$(ccgen) $(cflags) rk_6_3.c -o rk_6_3.pgo

rk_6_4.pgo: rk_6_4.c rk_6_4.h rk_6_3.h rk.h $(cfiles)
	$(ccgen) $(cflags) rk_6_4.c -o rk_6_4.pgo

ode.pgo: ode.c steps.h $(rkhfiles) $(cfiles)
	$(ccgen) $(cflags) ode.c -o ode.pgo

utils.gcda: ode-pgo script-tests.sh
	sh script-tests.sh

optimize.gcda: ode-pgo script-tests.sh
	sh script-tests.sh

rk.gcda: ode-pgo script-tests.sh
	sh script-tests.sh

rk_2_2.gcda: ode-pgo script-tests.sh
	sh script-tests.sh

rk_3_2.gcda: ode-pgo script-tests.sh
	sh script-tests.sh

rk_3_3.gcda: ode-pgo script-tests.sh
	sh script-tests.sh

rk_4_2.gcda: ode-pgo script-tests.sh
	sh script-tests.sh

rk_4_3.gcda: ode-pgo script-tests.sh
	sh script-tests.sh

rk_4_4.gcda: ode-pgo script-tests.sh
	sh script-tests.sh

rk_5_2.gcda: ode-pgo script-tests.sh
	sh script-tests.sh

rk_5_3.gcda: ode-pgo script-tests.sh
	sh script-tests.sh

rk_5_4.gcda: ode-pgo script-tests.sh
	sh script-tests.sh

rk_6_2.gcda: ode-pgo script-tests.sh
	sh script-tests.sh

rk_6_3.gcda: ode-pgo script-tests.sh
	sh script-tests.sh

rk_6_4.gcda: ode-pgo script-tests.sh
	sh script-tests.sh

ode.gcda: ode-pgo script-tests.sh
	sh script-tests.sh

ode.pdf: ode.tex Makefile
	pdflatex ode
	pdflatex ode
	pdflatex ode

clean:
	rm -rf ode *.{o,aux,toc,log} html latex

strip: ode
	strip ode
