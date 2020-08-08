./start.sh

for state in solid liquid gas
do

	echo
	echo ====================================
	echo
	echo 'STARTING main01.x for the '$state' state'

	#./start.sh without rm other states' output files
	cp config.fcc config.0
	cp seed.0 seed.in

	cp input.$state input.dat
	./main01.x $state

done
