for k in `seq 0 1`;
do
	for m in `seq 0 1`;
	do
		./ode-pgo tests/test-rk-2-2-$k-0-$m.xml
	done
done
for j in `seq 2 3`;
do
	for k in `seq 0 1`;
	do
		for l in `seq 0 1`;
		do
			for m in `seq 0 1`;
			do
				./ode-pgo tests/test-rk-3-$j-$k-$l-$m.xml
			done
		done
	done
done
for j in `seq 2 3`;
do
	for k in `seq 0 1`;
	do
		for l in `seq 0 1`;
		do
			for m in `seq 0 1`;
			do
				./ode-pgo tests/test-rk-4-$j-$k-$l-$m.xml
			done
		done
	done
done
for k in `seq 0 1`;
do
	for m in `seq 0 1`;
	do
		./ode-pgo tests/test-rk-4-4-$k-0-$m.xml
	done
done
for i in `seq 3 8`;
do
	k=`echo "$i-1" | bc`;
	for j in `seq 2 $k`;
	do
		echo "./ode-pgo tests/test-steps-$i-$j.xml"
		./ode-pgo tests/test-steps-$i-$j.xml
	done
done
for i in `seq 9 13`;
do
	for j in `seq 2 8`;
	do
		echo "./ode-pgo tests/test-steps-$i-$j.xml"
		./ode-pgo tests/test-steps-$i-$j.xml
	done
done
