for state in solid liquid gas
do

	echo
	echo ====================================
	echo
	echo 'STARTING main01.x for the '$state' state'

	./start.sh
	cp input.$state input.dat
	./main01.x $state

done
