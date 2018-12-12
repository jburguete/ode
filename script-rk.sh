# free=2
if ! test -s rk-2-2.mc; then
	./ode 1 2 2 64 10 64 0.5 0.5 7
fi
# free=5
if ! test -s rk-3-3.mc; then
	./ode 1 3 3 64 7 64 0.5 0.5 7
fi
# free=7
if ! test -s rk-3-2.mc; then
	./ode 1 3 2 64 6 64 0.5 0.5 7
fi
# free=8
if ! test -s rk-4-4.mc; then
	./ode 1 4 4 64 6 64 0.5 0.5 7
fi
# free=12
if ! test -s rk-4-3.mc; then
	./ode 1 4 3 64 4 64 0.5 0.5 7
fi
# free=14
if ! test -s rk-4-2.mc; then
	./ode 1 4 2 64 4 64 0.5 0.5 7
fi
