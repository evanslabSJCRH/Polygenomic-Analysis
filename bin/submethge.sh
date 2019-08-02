#!/bin/bash

cd ${3}/${4}
queue=${6}
conditional=""

while read probe
do
    cprobe=$(echo $probe | sed -e 's/\//_/g')
    
    conditional="${conditional}ended(${2}_${4}_B_4_${cprobe})&&"
    
    echo ${cprobe}
    
    bsub -P "${1}" -q ${queue}  -sp ${priority}  -J "${2}_${4}_D_1_${cprobe}[1-4]" -app R-2.14.0 -R "rusage[mem=30000]" -o log/${2}_${4}_\${LSB_JOBINDEX}_totxv_meth_ge.o R CMD BATCH --no-save --args --\${LSB_JOBINDEX} --${4} --$probe --TOTXV gemeth.R log/totxv_meth_ge_${4}_${cprobe}_\${LSB_JOBINDEX}.Rout
    bsub -P "${1}" -q ${queue}  -sp ${priority}  -J "${2}_${4}_D_2_${cprobe}[1-4]" -app R-2.14.0 -R "rusage[mem=30000]" -o log/${2}_${4}_\${LSB_JOBINDEX}_totxvi_meth_ge.o R CMD BATCH --no-save --args --\${LSB_JOBINDEX} --${4} --$probe --TOTXVI gemeth.R log/totxvi_meth_ge_${4}_${cprobe}_\${LSB_JOBINDEX}.Rout
 bsub -P "${1}" -q ${queue}  -sp ${priority}  -J "${2}_${4}_D_3_${cprobe}[1-4]" -app R-2.14.0 -R "rusage[mem=30000]" -o log/${2}_${4}_\${LSB_JOBINDEX}_all_meth_ge.o R CMD BATCH --no-save --args --\${LSB_JOBINDEX} --${4} --$probe --ALL gemeth.R log/all_meth_ge_${4}_${cprobe}_\${LSB_JOBINDEX}.Rout

    bsub -P "${1}" -q ${queue}  -sp ${priority}  -J "${2}_${4}_D_4_${cprobe}[1-4]" -app R-2.14.0 -R "rusage[mem=30000]" -o log/${2}_${4}_\${LSB_JOBINDEX}_meth_ge_meta.o -w "ended(${2}_${4}_D_1_${cprobe})&&ended(${2}_${4}_D_2_${cprobe})" R CMD BATCH --no-save --args --\${LSB_JOBINDEX} --${4} --$probe gemethmeta.R log/meth_ge_meta_${4}_${cprobe}_\${LSB_JOBINDEX}.Rout
    
    bsub -P "${1}" -q ${queue}  -sp ${priority}  -J "${2}_${4}_D_5_${cprobe}[1-4]" -app R-2.14.0 -R "rusage[mem=30000]" -o log/${2}_${4}_\${LSB_JOBINDEX}_meth_ge_meta_sub.o -w "ended(${2}_${4}_D_3_${cprobe})" R CMD BATCH --no-save --args --\${LSB_JOBINDEX} --${4} --$probe gemethmetasub.R log/meth_ge_meta_sub_${4}_${cprobe}_\${LSB_JOBINDEX}.Rout


done < meth_sig_probes.tsv

#echo ${conditional}

conditional="${conditional%??}"

bsub -P "${1}" -q ${queue}  -sp ${priority}  -J "${2}_${4}_D_6[1]" -R "rusage[mem=30000]" -o log/${2}_${4}_\${LSB_JOBINDEX}_meth_ge_meta_sub_result_cat.o -w "${conditional}" ${5}/checkconcat.sh ${3}/${4}/ge_meth_meta_sub ${3}/${4}/ge_meth_meta_sub.tsv
