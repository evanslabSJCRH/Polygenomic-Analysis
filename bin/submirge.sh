#!/bin/bash
priority=${6}
cd ${3}/${4}

conditional=""
queue=${6}

while read probe
do
    cprobe=$(echo $probe | sed -e 's/\//_/g')
    
    conditional="${conditional}ended(${2}_${4}_C_3_${cprobe})&&"
    
    echo ${cprobe}
    
    bsub -P "${1}" -q ${queue}  -sp ${priority} -J "${2}_${4}_C_1_${cprobe}" \
        -app R-2.14.0 -R "rusage[mem=2000]" \
        -o log/${2}_${4}_${cprobe}_totxv_mir_ge.o \
        R CMD BATCH --no-save --args --$probe --${4} \
        --TOTXV gemir.R \
        log/totxv_mir_ge_${4}_${cprobe}.Rout

    bsub -P "${1}" -q ${queue}  -sp ${priority} -J "${2}_${4}_C_2_${cprobe}" \
        -app R-2.14.0 -R "rusage[mem=2000]" \
        -o log/${2}_${4}_${cprobe}_totxvi_mir_ge.o \
        R CMD BATCH --no-save --args --$probe --${4} \
        --TOTXVI gemir.R \
        log/totxvi_mir_ge_${4}_${cprobe}.Rout

    bsub -P "${1}"  -q ${queue}  -sp ${priority} -J "${2}_${4}_C_2all_${cprobe}" \
        -app R-2.14.0 -R "rusage[mem=2000]" \
        -o log/${2}_${4}_${cprobe}_all_mir_ge.o \
        R CMD BATCH --no-save --args --$probe --${4} \
        --ALL gemir.R \
        log/all_mir_ge_${4}_${cprobe}.Rout

    bsub -P "${1}"  -q ${queue}  -sp ${priority} -J "${2}_${4}_C_3_${cprobe}" \
        -app R-2.14.0 -R "rusage[mem=30000]" \
        -o log/${2}_${4}_${cprobe}_mir_ge_meta.o \
        -w "ended(${2}_${4}_C_1_${cprobe})&&ended(${2}_${4}_C_2_${cprobe})&&ended(${2}_${4}_C_2all_${cprobe})" \
        R CMD BATCH --no-save --args --\${LSB_JOBINDEX} --${4} \
        --$probe gemirmeta.R \
        log/mir_ge_meta_${4}_${cprobe}.Rout
    
done < mir_sig_probes.tsv

conditional="${conditional%??}"
echo ${conditional}

bsub -P "${1}" -q ${queue}  -sp ${priority} -J "${2}_${4}_C_5" -R "rusage[mem=30000]" \
    -o log/${2}_${4}_ge_mir_meta_sub_result_cat.o \
    -w "${conditional}" \
    ${5}/checkconcat.sh ${3}/${4}/ge_mir_meta ${3}/${4}/ge_mir_meta.tsv
