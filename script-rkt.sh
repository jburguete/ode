# free=0
if ! test -s rk-2-2.mc; then
	./odet 3 2 2 1 1 1 1 1 7
fi
# free=1
if ! test -s rk-3-3.mc; then
	./odet 3 3 3 64 1000000 64 0.5 0.5 7
fi
# free=3
if ! test -s rk-3-2.mc; then
	./odet 3 3 2 64 100 64 0.5 0.5 7
fi
# free=1
if ! test -s rk-4-4.mc; then
	./odet 3 4 4 64 1000000 64 0.5 0.5 7
fi
# free=5
if ! test -s rk-4-3.mc; then
	./odet 3 4 3 64 16 64 0.5 0.5 7
fi
# free=7
if ! test -s rk-4-2.mc; then
	./odet 3 4 2 64 4 64 0.5 0.5 7
fi
# free=6
if ! test -s rk-5-4.mc; then
	./odet 3 5 4 64 10 64 0.5 0.5 7
fi
# free=10
if ! test -s rk-5-3.mc; then
	./odet 3 5 3 64 4 64 0.5 0.5 7
fi
# free=12
if ! test -s rk-5-2.mc; then
	./odet 3 5 2 64 3 64 0.5 0.5 7
fi
