
module load EBModules
module load kallisto/0.46.1-foss-2019b

ulimit -n 10000

kallisto quant -i $2 -t 4 \
-o ${1} \
-b 100 $3 $4


