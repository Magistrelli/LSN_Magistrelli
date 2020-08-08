./start.sh

for state in solid liquid gas
do
	echo
	echo ====================================
	echo
	echo 'STARTING main03.x for the '$state' state'

	#./start.sh without rm other states' output files	
	rm -rf old.0
	cp config.fcc config.0
	cp seed.0 seed.in

	cp inputNVE.$state input.dat
	./main03.x $state 10
done

for i in 20 30
do
	echo
	echo ====================================
	echo
	echo 'STARTING main03.x for the gas state with more ('$i') restarts'
	
	rm -rf old.0
	cp config.fcc config.0
	cp seed.0 seed.in

	cp inputNVE.gas input.dat
	./main03.x gas $i
done
