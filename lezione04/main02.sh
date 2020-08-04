./start.sh

for state in solid liquid gas
do

	echo
	echo ====================================
	echo
	echo 'STARTING main01.x for the '$state' state'

	#./start.sh without rm other states' output files	
	rm -rf old.0
	cp config.fcc config.0
	cp seed.0 seed.in

	#executing main01.x for the first time
	cp input_short.$state input.dat
	./main01.x $state 0

	#restarting other 7 time main01.x
	for i in {1..7}
	do
		echo
		echo ====================================
		echo
		echo 'Restarting '$state' for the '$i'-th time'
		./restart.sh
		./main01.x $state $i
	done

	echo
	echo ====================================
	echo
	echo 'STARTING main02.x for the '$state' state'
	./restart.sh
	cp input.$state input.dat
	./main02.x $state
done
