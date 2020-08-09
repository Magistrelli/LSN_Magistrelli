for state in solid liquid gas
do
	echo
	echo ====================================
	echo
	echo 'STARTING main03.x for the '$state' state'
	
	rm -rf old.0
	./start.sh

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
	./start.sh

	cp inputNVE.gas input.dat
	./main03.x gas $i
done
