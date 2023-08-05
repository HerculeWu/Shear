if [ $# -le 1 ]
then
    echo "Usage:"
    echo "sh ./run.sh [number for mpiexec] [number of runs for each epsilonh] [delta] [A]"
else
    for epsilonh in 0.1 0.2; do
        run=0
        while [ $run -le $2 ]
        do
            echo "epsilonh: $epsilonh run: $run"
            python3 ./generator.py 25 30 ./data.$epsilonh.$run.in
            mpiexec -n $1 ./shear ./data.$epsilonh.$run.in $3 $4 $epsilonh
            run=$(( $run+1 ))
            # gether
            touch out.$epsilonh.$run
            rank=0
            while [ $rank -le $(( $1-1 )) ]
            do
                cat ./out_$rank.txt >> out.$epsilonh.$run
                rank=$(( $rank+1 ))
            done
        done
    done
fi