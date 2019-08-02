#!/bin/bash
priority="100"
queue=${6}


cd ${3}/${4}
#TODO(spaugh):Add job that sleeps for calculated time \
#of submission and make all A_1 jobs depend on that job\
#so that all jobs are pending when final concat job is submitted.

conditional=""

while read probe
do
    cprobe=$(echo $probe | sed -e 's/\//_/g')
    
    conditional="${conditional}ended(${2}_${4}_A_4_${cprobe})&&"
    
    echo ${cprobe}
    
    bsub -P "${1}" -q ${queue} -sp ${priority} -J "${2}_${4}_A_1_${cprobe}[1-4]" \
        -app R-2.14.0 -R "rusage[mem=30000]" \
        -oo "log/${2}_${4}_${cprobe}_%I_totxv_snp_ge.o" \
        R CMD BATCH --no-save --args --\${LSB_JOBINDEX} \
        --${4} --$probe --TOTXV gesnp.R \
        log/totxv_snp_ge_${4}_${cprobe}_\${LSB_JOBINDEX}.Rout
    
    bsub -P "${1}" -q ${queue}  -sp ${priority} -J "${2}_${4}_A_2_${cprobe}[1-4]" \
        -app R-2.14.0 -R "rusage[mem=30000]" \
        -oo "log/${2}_${4}_${cprobe}_%I_totxvi_snp_ge.o" \
        R CMD BATCH --no-save --args --\${LSB_JOBINDEX} \
        --${4} --$probe --TOTXVI gesnp.R \
        log/totxvi_snp_ge_${4}_${cprobe}_\${LSB_JOBINDEX}.Rout

  bsub -P "${1}" -q ${queue} -sp ${priority} -J "${2}_${4}_A_3_${cprobe}[1-4]" \
        -app R-2.14.0 -R "rusage[mem=30000]" \
        -oo "log/${2}_${4}_${cprobe}_%I_all_snp_ge.o" \
        R CMD BATCH --no-save --args --\${LSB_JOBINDEX} \
        --${4} --$probe --ALL gesnp.R \
        log/all_snp_ge_${4}_${cprobe}_\${LSB_JOBINDEX}.Rout

    bsub -P "${1}" -q ${queue} -sp ${priority} -J "${2}_${4}_A_4_${cprobe}[1-4]" \
        -app R-2.14.0 -R "rusage[mem=30000]" \
        -oo "log/${2}_${4}_${cprobe}_%I_snp_ge_meta.o" \
        -w "ended(${2}_${4}_A_1_${cprobe})&&ended(${2}_${4}_A_2_${cprobe})&&ended(${2}_${4}_A_3_${cprobe})" \
        R CMD BATCH --no-save --args --\${LSB_JOBINDEX} \
        --${4} --$probe gesnpmeta.R \
        log/snp_ge_meta_${4}_${cprobe}_\${LSB_JOBINDEX}.Rout

    bsub -P "${1}" -q ${queue}  -sp ${priority} -J "${2}_${4}_A_5_${cprobe}[1-4]" \
        -app R-2.14.0 -R "rusage[mem=30000]" \
        -oo "log/${2}_${4}_${cprobe}_%I_snp_ge_meta_sub.o" \
        -w "ended(${2}_${4}_A_4_${cprobe})" \
        R CMD BATCH --no-save --args --\${LSB_JOBINDEX} \
        --${4} --$probe gesnpmetasub.R \
        log/snp_ge_meta_sub_${4}_${cprobe}_\${LSB_JOBINDEX}.Rout

done < ge_snp_sig_probes.tsv

conditional="${conditional%??}"

echo ${conditional}

bsub -P "${1}" -q ${queue}  -sp ${priority} -J "${2}_${4}_A_6[1]" \
    -R "rusage[mem=30000]" -oo "log/${2}_${4}_%I_snp_ge_meta_sub_result_cat.o" \
    -w "${conditional}" ${5}/concat.sh \
    ${3}/${4}/ge_snp_meta_sub ${3}/${4}/ge_snp_meta_sub.tsv


#bsub -P $project -sp ${priority} -J ${timestamp}_${drug}_ge_snp_cutoff[1] -app R-2.14.0 -R "rusage[mem=30000]" R CMD BATCH --no-save --args ${drug}_ge_snp_cutoff.R log/${drug}_ge_snp_cutoff\$LSB_JOBINDEX.Rout	
