for j in 1 2
do
	cp seed.out seed.in
	for i in {0..20}
	do
		cp finalConf/config.metro.$i.final config.0
		./main01.x $i $j
	done
done
