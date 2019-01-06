./ode-pgo tests/test-rk-2-2-0-0-0.xml
./ode-pgo tests/test-rk-2-2-0-0-1.xml
./ode-pgo tests/test-rk-2-2-1-0-0.xml
./ode-pgo tests/test-rk-2-2-1-0-1.xml
./ode-pgo tests/test-rk-3-2-0-0-0.xml
./ode-pgo tests/test-rk-3-2-0-0-1.xml
./ode-pgo tests/test-rk-3-2-1-0-0.xml
./ode-pgo tests/test-rk-3-2-1-0-1.xml
for i in `seq 3 8`;
do
	k=`echo "$i-1" | bc`;
	echo $k
	for j in `seq 2 $k`;
	do
		echo "./ode-pgo tests/test-steps-$i-$j.xml"
		./ode-pgo tests/test-steps-$i-$j.xml
	done
done
for i in `seq 9 12`;
do
	for j in `seq 2 8`;
	do
		echo "./ode-pgo tests/test-steps-$i-$j.xml"
		./ode-pgo tests/test-steps-$i-$j.xml
	done
done
