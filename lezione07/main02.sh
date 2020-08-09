for state in solid liquid gas
do

	echo
	echo ====================================
	echo
	echo 'STARTING main02.x for the '$state' state'

	./start.sh
	cp input.$state input.dat
	./main02.x $state 0

done
