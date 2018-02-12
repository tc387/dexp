#!/bin/bash/

for rep in 1 ;
do
    let xseed=12341-4*$rep  # random seed

    for nlig in 256
    do
	for rcol in 50
	do
	    for maxnbonds in 24 36
	    do
		let maxbonddist=30
		
		for form in m t r # mobile, annealed, thomson, random
		do
		    ndir='nsim_dexp_f'$form'_mnb'$maxnbonds'_nl'$nlig'_rc'$rcol'_r'$rep
		    
		    if [[ "$form" == "m" ]]
		    then
			moblig='.true.'
			nanneal=0
			readinitconf='.false.'
		    elif [[  "$form" == "a" ]]
		    then
			moblig='.false.'
			nanneal=100000
			readinitconf='.false.'
		    elif [[  "$form" == "t" ]]
		    then
			moblig='.false.'
			nanneal=0
			readinitconf='.true.'
		    elif [[  "$form" == "r" ]]
		    then
			moblig='.false.'
			nanneal=0
			readinitconf='.false.'
		    fi
		    
		    echo $moblig
		    echo $nanneal
		    echo $readinitconf
		    
		    cat input-wlb-master.dat | sed -e 's/XSEED11/'$xseed'/' |  sed -e 's/TNCH11/'$nlig'/g' | sed -e 's/RCOL11/'$rcol'/' | sed -e 's/MAXBD11/'$maxbonddist'/'  | sed -e 's/MOBLIG11/'$moblig'/'  | sed -e 's/NANNEAL11/'$nanneal'/'  | sed -e 's/READINITC11/'$readinitconf'/' | sed -e 's/MAXNB11/'$maxnbonds'/' > input-parwlb.dat
                    # cat qsub-master.sh | sed -e 's/JOBNAME11/limem_eps'$eps'_th33'$th33'/' > qsub_local.sh
		    mkdir  $ndir
		    cp input-parwlb.dat a.out thomson$nlig'.xyz' $ndir
		    cd $ndir
		    ./a.out > screen.txt &
		    cd ..
		    
		done
	    done
	done
    done
done
		    
