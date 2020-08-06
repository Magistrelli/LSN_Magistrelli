./clean.sh

for met in metro gibbs
do
	for h in input input_h
	do
		rm -rf config.0
		cp seed.0 seed.in
		cp $h.$met input.dat

		for i in {0..20}
		do
			./main01.x $i 0
		done
	done
done
