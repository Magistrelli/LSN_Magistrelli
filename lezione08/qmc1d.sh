#PIGS
for tau in 05 1 2 3
do
	cp input/input_t$tau.pigs input.dat
	./qmc1d.x	#variational trial wave function

	for res in probability potential kinetic
	do
		cp $res.out results/$res.var_t$tau.out
	done;echo;echo	

done


for tau in 3 5 6
do
	cp input/input_t$tau.pigs input.dat
	./qmc1d_const.x	#constant trial wave function

	for res in probability potential kinetic
	do
		cp $res.out results/$res.const_t$tau.out
	done;echo;echo

done


#PIMC
for temp in 0_25 1_25 5 50
do
	cp input/input_$temp.pimc input.dat
	./qmc1d.x

	for res in probability potential kinetic
	do
		cp $res.out results/$res.T$temp.out
	done;echo;echo

done

for res in probability potential kinetic
do
	rm -rf $res.out
done
