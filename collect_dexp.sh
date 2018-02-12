#!/bin/bash/

deleps=0.1

for nlig in 256
do
    for rcol in 30 40 50
    do
	for maxnbonds in 24
	do

	    for form in m a t r # mobile, annealed, thomson, random
	    do
		for rep in 1 ;
		do
		    ndir='nsim_dexp_f'$form'_mnb'$maxnbonds'_nl'$nlig'_rc'$rcol'_r'$rep
		    pyoutfile='../out_Fofeps_f'$form'_mnb'$maxnbonds'_nl'$nlig'_rc'$rcol'_r'$rep'.dat'

		    cd $ndir
		    for ifile in $(ls | grep test_WLFE ) ;
		    do
			lastoutfile=$ifile
		    done
		    echo "" 
		    echo $ndir
		    echo $lastoutfile
		    python3 ../get_Fofeps_dexp.py $lastoutfile $pyoutfile $deleps
		    
	
		    cd ..
		done
	    done
	done
    done
done
		    
