#!/bin/bash

cd ${3}/${4}

conditional=""

while read probe
do
    cprobe=$(echo $probe | sed -e 's/\//_/g')

    conditional="${conditional}ended(${2}_${4}_B_4_${cprobe})&&"

    echo ${cprobe}

    bsub -P "${1}" -q ${6} -J "${2}_${4}_B_1_${cprobe}[1-4]" \
        -app R-2.14.0 -R "rusage[mem=30000]" \
        -o log/${2}_${4}_${LSB_JOBINDEX}_totxv_cn_ge.o \
        R CMD BATCH --no-save --args --\${LSB_JOBINDEX} \
        --${4} --$probe --TOTXV gecn.R \
        log/totxv_cn_ge_${4}_${cprobe}_\${LSB_JOBINDEX}.Rout

    bsub -P "${1}" -q ${6} -J "${2}_${4}_B_2_${cprobe}[1-4]" \
        -app R-2.14.0 -R "rusage[mem=30000]" \
        -o log/${2}_${4}_${LSB_JOBINDEX}_totxvi_cn_ge.o \
        R CMD BATCH --no-save --args --\${LSB_JOBINDEX} \
        --${4} --$probe --TOTXVI gecn.R \
        log/totxvi_cn_ge_${4}_${cprobe}_\${LSB_JOBINDEX}.Rout

    bsub -P "${1}" -q ${6} -J "${2}_${4}_B_3_${cprobe}[1-4]" \
        -app R-2.14.0 -R "rusage[mem=30000]" \
        -o log/${2}_${4}_${LSB_JOBINDEX}_all_cn_ge.o \
        R CMD BATCH --no-save --args --\${LSB_JOBINDEX} \
        --${4} --$probe --ALL gecn.R \
        log/all_cn_ge_${4}_${cprobe}_\${LSB_JOBINDEX}.Rout

    bsub -P "${1}" -q ${6} -J "${2}_${4}_B_4_${cprobe}[1-4]" \
        -app R-2.14.0 -R "rusage[mem=30000]" \
        -o log/${2}_${4}_${LSB_JOBINDEX}_cn_ge_meta.o \
        -w "ended(${2}_${4}_B_1_${cprobe})&&ended(${2}_${4}_B_2_${cprobe})&&ended(${2}_${4}_B_3_${cprobe})" \
        R CMD BATCH --no-save --args --\${LSB_JOBINDEX} --${4} \
        --$probe gecnmeta.R \
        log/cn_ge_meta_${4}_${cprobe}_\${LSB_JOBINDEX}.Rout

    bsub -P "${1}" -q ${6} -J "${2}_${4}_B_5_${cprobe}[1-4]" \
        -app R-2.14.0 -R "rusage[mem=30000]" \
        -o log/${2}_${4}_${LSB_JOBINDEX}_cn_ge_meta_sub.o \
        -w "ended(${2}_${4}_B_4_${cprobe})" \
        R CMD BATCH --no-save --args --\${LSB_JOBINDEX} --${4} \
        --$probe gecnmetasub.R log/cn_ge_meta_sub_${4}_${cprobe}_\${LSB_JOBINDEX}.Rout

done < ge_sig_probes_cn.tsv

echo ${conditional}

conditional="${conditional%??}"

bsub -P "${1}" -q ${6} -J "${2}_${4}_B_6[1]" \
    -R "rusage[mem=30000]" -o log/${2}_${4}_${LSB_JOBINDEX}_cn_ge_meta_sub_result_cat.o \
    -w "${conditional}" \
    ${5}/checkconcat.sh ${3}/${4}/ge_cn_meta_sub ${3}/${4}/ge_cn_meta_sub.tsv
