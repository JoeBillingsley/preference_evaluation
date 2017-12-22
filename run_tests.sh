#!/bin/bash

command > tests.out 2>test.error

for dim_idx in {0..4}
do
	for run_idx in {0..20}
	do
		for test_idx in {0..7}  
		do  	
			for algo_idx in {0..7}
			do
				if [ $test_idx -lt 4 ] && [ $dim_idx -gt 0 ]
				then
					break
				fi

				./Samaritan $test_idx $algo_idx $dim_idx $run_idx
			done
		done
	done
done
