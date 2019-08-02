#!/bin/bash
#
#Runs the drworkflow pipeline.
################################################################################

set -e

source config.rc

timestamp=`date +%Y%m%d%H%M%S%N`
currdate=`date +%Y%m%d_%H%M%S`

echo $wd

echo $timestamp

# Prestage the files
datadir="/drworkflow_data"
ln -sf ${datadir}/20111103_SOM_GERM_SNP6_TOTXVI.csv \
    ./bin/20111103_SOM_GERM_SNP6_TOTXVI.csv
ln -sf ${datadir}/all.CNState.RData ./bin/all.CNState.RData
ln -sf ${datadir}/all.CNState.smooth.RData ./bin/all.CNState.smooth.RData
ln -sf ${datadir}/CNState.RData ./bin/CNState.RData 
ln -sf ${datadir}/dutch.RData ./bin/dutch.RData
ln -sf ${datadir}/human_predictions_S_0_aug2010.txt \
    ./bin/human_predictions_S_0_aug2010.txt
ln -sf ${datadir}/human_predictions_S_C_aug2010.txt \
    ./bin/human_predictions_S_C_aug2010.txt
ln -sf ${datadir}/SmoothSignal.RData ./bin/SmoothSignal.RData
ln -sf ${datadir}/2013-11-06.dutch.dxbm.ii.mas5.probe.log2.RData ./bin/2013-11-06.dutch.dxbm.ii.mas5.probe.log2.RData

for drug in "${drugs[@]}"
do
    echo $drug
    
    mkdir -p "$wd/$drug"
    mkdir -p "$wd/$drug/log"
    cd "$wd/$drug"

################################################################################
    # Stage the files in each drug output directory

    
    
    mkdir -p {totxv_snp_lc50,totxvi_snp_lc50,all_snp_lc50,snp_lc50_meta}
    mkdir -p {totxv_cn_lc50,totxvi_cn_lc50,all_cn_lc50,cn_lc50_meta}
    
    cp ${execbin}/cnlc50.R ${wd}/${drug}/cnlc50.R
    cp ${execbin}/genecnlc50.R ${wd}/${drug}/genecnlc50.R
    cp ${execbin}/genecncutoff.R ${wd}/${drug}/genecncutoff.R
    cp ${execbin}/drugcatadjust.R ${wd}/${drug}/drugcatadjust.R
    cp ${execbin}/U133_gene_pos.txt ${wd}/${drug}/U133_gene_pos.txt
    cp ${execbin}/mir_all_small.csv ${wd}/${drug}/mir_all_small.csv

    ln -sf ${execbin}/SmoothSignal.RData ${wd}/${drug}/SmoothSignal.RData
    ln -sf ${execbin}/all.CNState.RData ${wd}/${drug}/all.CNState.RData
    ln -sf ${execbin}/CNState.RData ${wd}/${drug}/CNState.RData
    ln -sf ${execbin}/all.CNState.smooth.RData \
        ${wd}/${drug}/all.CNState.smooth.RData
    cp -u ${execbin}/20111103_SOM_GERM_SNP6_TOTXVI.csv \
        $wd/$drug/20111103_SOM_GERM_SNP6_TOTXVI.csv
    cp ${execbin}/cnlc50merge.R $wd/$drug/cnlc50merge.R
    cp ${execbin}/cnlc50meta.R $wd/$drug/cnlc50meta.R


    cp ${execbin}/Ip_Select-Optimal-Pval-Threshold_7-22-12.R \
        ${wd}/${drug}/Ip_Select-Optimal-Pval-Threshold_7-22-12.R    

    ln -sf ${execbin}/human_predictions_S_C_aug2010.txt \
        ${wd}/${drug}/human_predictions_S_C_aug2010.txt 
    ln -sf ${execbin}/human_predictions_S_0_aug2010.txt \
        ${wd}/${drug}/human_predictions_S_0_aug2010.txt
    ln -sf ${execbin}/miranno.csv $wd/$drug/miranno.csv
    ln -sf ${execbin}/miranno_new.csv $wd/$drug/miranno_new.csv
    ln -sf ${execbin}/mir_all.csv $wd/$drug/mir_all.csv
    cp -u ${execbin}/Sweave.sty ${wd}/${drug}/Sweave.sty

    ln -sf ${execbin}/dutch.RData $wd/$drug/dutch.RData
    ln -sf ${execbin}/nl_profile051405.txt $wd/$drug/nl_profile051405.txt

    create_sweave_dispatch()
    {
        cp ${execbin}/${1}.Rnw $wd/$drug/${1}.Rnw
        echo -e "r.lib <- '/home/rautry/drworkflow_Rlib'" > $wd/$drug/${1}.R
        echo -e "require (filehash, lib.loc=r.lib)" >> $wd/$drug/${1}.R
        echo -e "require (digest, lib.loc=r.lib)" >> $wd/$drug/${1}.R
        echo -e "require (stashR, lib.loc=r.lib)" >> $wd/$drug/${1}.R
        echo -e "require(cacheSweave, lib.loc=r.lib)" >> $wd/$drug/${1}.R
        echo -e "setCacheDir('./cache${1}/')" >> $wd/$drug/${1}.R
        echo -e "Sweave ('${1}.Rnw', driver=cacheSweaveDriver)" >> \
            $wd/$drug/${1}.R
    }

    create_sweave_dispatch gelc50
    create_sweave_dispatch mirlc50
    create_sweave_dispatch methlc50
    

#################################################################################

if [ "${drug}" == "6TG" ]; then

    cp ${execbin}/6TG_LC50_AllCohorts_022414.csv "${wd}/${drug}/6TG_LC50_AllCohorts_022414.csv"
    ln -sf ${execbin}/2013-11-06.dutch.dxbm.ii.mas5.probe.log2.RData ${wd}/${drug}/2013-11-06.dutch.dxbm.ii.mas5.probe.log2.RData

fi
 
#################################################################################

if [ "${drug}" == "6MP" ]; then

    cp ${execbin}/6MP_LC50_AllCohorts_022414.csv "${wd}/${drug}/6MP_LC50_AllCohorts_022414.csv"
    ln -sf ${execbin}/2013-11-06.dutch.dxbm.ii.mas5.probe.log2.RData ${wd}/${drug}/2013-11-06.dutch.dxbm.ii.mas5.probe.log2.RData
fi

#################################################################################
if [ "${runsnp}" == "TRUE" ]; then

    cp ${execbin}/snpcutoff.R "${wd}/${drug}/${drug}_snpcutoff.R"
    cp ${execbin}/snplc50.R "$wd/$drug/snplc50.R"
    cp ${execbin}/snplc50merge.R "$wd/$drug/snplc50merge.R"
    cp ${execbin}/snplc50meta.R "$wd/$drug/snplc50meta.R"
    mkdir -p {imputedxv_snp_lc50,imputedxvi_snp_lc50}

 
    bsub -P $project -q ${queue} -sp ${priority} -J ${timestamp}_${drug}_01_01_01[1-200] \
        -oo "log/${timestamp}_${drug}_01_01_01_%I.o" -R "rusage[mem=30000]" -M 30000 \
        -app R-2.14.0 R CMD BATCH --no-save --args --\$LSB_JOBINDEX --$drug \
        --TOTXV snplc50.R log/${timestamp}_${drug}_01_01_01_\$LSB_JOBINDEX.Rout

    bsub -P $project -q ${queue} -sp ${priority} -J ${timestamp}_${drug}_01_01_02[1] \
        -oo "log/${timestamp}_${drug}_01_01_02.o" -R "rusage[mem=30000]" -M 30000 \
        -w "ended(${timestamp}_${drug}_01_01_01)" $execbin/checkconcat.sh \
        $wd/$drug/totxv_snp_lc50 $wd/$drug/totxv_snp_lc50.tsv

    bsub -P $project -q ${queue} -sp ${priority} -J ${timestamp}_${drug}_01_01_03[1-200] \
        -oo "log/${timestamp}_${drug}_01_01_03_%I.o" -R "rusage[mem=30000]" -M 30000 \
        -app R-2.14.0 R CMD BATCH --no-save --args --\$LSB_JOBINDEX --$drug \
        --TOTXVI snplc50.R log/${timestamp}_${drug}_01_01_03_\$LSB_JOBINDEX.Rout

    bsub -P $project -q ${queue} -sp ${priority} -J ${timestamp}_${drug}_01_01_04[1] \
        -oo "log/${timestamp}_${drug}_01_01_04.o" \
        -w "ended(${timestamp}_${drug}_01_01_03)" -R "rusage[mem=30000]" -M 30000 \
        $execbin/checkconcat.sh $wd/$drug/totxvi_snp_lc50 $wd/$drug/totxvi_snp_lc50.tsv


    bsub -P $project -q ${queue} -sp ${priority} -J ${timestamp}_${drug}_01_01_05[1-200] \
        -app R-2.14.0 -oo "log/${timestamp}_${drug}_01_01_05_%I.o" \
        -R "rusage[mem=30000]" -M 30000 R CMD BATCH --no-save --args --\$LSB_JOBINDEX \
        --$drug --ALL snplc50.R log/${timestamp}_${drug}_01_01_05_\$LSB_JOBINDEX.Rout

    bsub -P $project -q ${queue} -sp ${priority} -J ${timestamp}_${drug}_01_01_06[1] \
        -oo "log/${timestamp}_${drug}_01_01_06.o" \
        -w "ended(${timestamp}_${drug}_01_01_05)" -R "rusage[mem=30000]" -M 30000 \
        $execbin/checkconcat.sh $wd/$drug/all_snp_lc50 $wd/$drug/all_snp_lc50.tsv


    bsub -P $project -q ${queue} -sp ${priority} -J ${timestamp}_${drug}_01_01_07[1-285,288,289] \
        -oo "log/${timestamp}_${drug}_01_01_07_%I.o" -R "rusage[mem=30000]" -M 30000 \
        -app R-2.14.0 R CMD BATCH --no-save --args --\$LSB_JOBINDEX --$drug \
        --imputedxv snplc50.R log/${timestamp}_${drug}_01_01_07_\$LSB_JOBINDEX.Rout

    bsub -P $project -q ${queue} -sp ${priority} -J ${timestamp}_${drug}_01_01_08[1] \
        -oo "log/${timestamp}_${drug}_01_01_08.o" -R "rusage[mem=30000]" -M 30000 \
        -w "ended(${timestamp}_${drug}_01_01_07)" $execbin/checkconcat.sh \
        $wd/$drug/imputedxv_snp_lc50 $wd/$drug/imputedxv_snp_lc50.tsv

    bsub -P $project -q ${queue} -sp ${priority} -J ${timestamp}_${drug}_01_01_09[1-289] \
        -oo "log/${timestamp}_${drug}_01_01_09_%I.o" -R "rusage[mem=30000]" -M 30000 \
        -app R-2.14.0 R CMD BATCH --no-save --args --\$LSB_JOBINDEX --$drug \
        --imputedxvi snplc50.R log/${timestamp}_${drug}_01_01_09_\$LSB_JOBINDEX.Rout

    bsub -P $project -q ${queue} -sp ${priority} -J ${timestamp}_${drug}_01_01_10[1] \
        -oo "log/${timestamp}_${drug}_01_01_10.o" -R "rusage[mem=30000]" -M 30000 \
        -w "ended(${timestamp}_${drug}_01_01_09)" $execbin/checkconcat.sh \
        $wd/$drug/imputedxvi_snp_lc50 $wd/$drug/imputedxvi_snp_lc50.tsv

    bsub -P $project -q ${queue} -sp ${priority} -J ${timestamp}_${drug}_01_01_11[1] \
        -app R-2.14.0 -oo "log/${timestamp}_${drug}_01_01_11.o" \
        -w "ended(${timestamp}_${drug}_01_01_02)&&ended(${timestamp}_${drug}_01_01_04)&&ended(${timestamp}_${drug}_01_01_06)" \
        -R "rusage[mem=2000]" -M 2000 R CMD BATCH --no-save --args --$drug snplc50merge.R \
        log/${timestamp}_${drug}_01_01_11_\$LSB_JOBINDEX.Rout

    bsub -P $project -q ${queue} -sp ${priority} -J ${timestamp}_${drug}_01_01_12[1-200] \
        -app R-2.14.0 -oo "log/${timestamp}_${drug}_01_01_12_%I.o" \
        -w "ended(${timestamp}_${drug}_01_01_11)" -R "rusage[mem=2000]" -M 2000 \
        R CMD BATCH --no-save --args --\$LSB_JOBINDEX --$drug snplc50meta.R \
        log/${timestamp}_${drug}_01_01_12_\$LSB_JOBINDEX.Rout

    bsub -P $project -q ${queue} -sp ${priority} -J ${timestamp}_${drug}_01_01_13[1] \
        -oo "log/${timestamp}_${drug}_01_01_13.o" \
        -R "rusage[mem=30000]" -M 30000 -w "ended(${timestamp}_${drug}_01_01_12)" \
        $execbin/checkconcat.sh $wd/$drug/snp_lc50_meta $wd/$drug/snp_lc50.tsv

    bsub -P $project -q ${queue} -sp ${priority} -J ${timestamp}_${drug}_01_01_14[1] \
        -oo "log/${timestamp}_${drug}_01_01_14.o" \
        -app R-2.14.0 -R "rusage[mem=30000]" -M 30000 \
        -w "ended(${timestamp}_${drug}_01_01_13)" \
        R CMD BATCH --no-save --args ${drug}_snpcutoff.R \
        log/${timestamp}_${drug}_01_01_14.Rout

fi
#################################################################################    



#################################################################################    

if [[ "${runcn}" == "TRUE" ]]; then

  # cp ${execbin}/cncutoff.R ${wd}/${drug}/${drug}_cncutoff.R
   cp ${execbin}/genecncutoff.R ${wd}/${drug}/${drug}_cncutoff.R

    bsub -P $project -q ${queue} -sp ${priority} -J ${timestamp}_${drug}_01_02_01[1-200] \
        -app R-2.14.0 -oo "log/${timestamp}_${drug}_01_02_01_%I.o" \
        -R "rusage[mem=30000]" -M 30000 R CMD BATCH --no-save --args --\$LSB_JOBINDEX \
        --$drug --TOTXV genecnlc50.R log/${timestamp}_${drug}_01_02_01_\$LSB_JOBINDEX.Rout

    bsub -P $project -q ${queue} -sp ${priority} -J ${timestamp}_${drug}_01_02_02[1] \
        -oo "log/${timestamp}_${drug}_01_02_02.o" \
        -R "rusage[mem=30000]" -M 30000 -w "ended(${timestamp}_${drug}_01_02_01)" \
        $execbin/checkconcat.sh $wd/$drug/totxv_cn_lc50 $wd/$drug/totxv_cn_lc50.tsv


    bsub -P $project -q ${queue} -sp ${priority} -J ${timestamp}_${drug}_01_02_03[1-200] \
        -app R-2.14.0 -oo "log/${timestamp}_${drug}_01_02_03_%I.o" \
        -R "rusage[mem=30000]" -M 30000 R CMD BATCH --no-save --args --\$LSB_JOBINDEX \
        --$drug --TOTXVI genecnlc50.R log/${timestamp}_${drug}_01_02_03_\$LSB_JOBINDEX.Rout

    bsub -P $project -q ${queue} -sp ${priority} -J ${timestamp}_${drug}_01_02_04[1] \
        -oo "log/${timestamp}_${drug}_01_02_04.o" \
        -R "rusage[mem=30000]" -M 30000 -w "ended(${timestamp}_${drug}_01_02_03)" \
        $execbin/checkconcat.sh $wd/$drug/totxvi_cn_lc50 $wd/$drug/totxvi_cn_lc50.tsv


    bsub -P $project -q ${queue} -sp ${priority} -J ${timestamp}_${drug}_01_02_05[1-200] \
        -app R-2.14.0 -oo "log/${timestamp}_${drug}_01_02_05_%I.o" \
        -R "rusage[mem=30000]" -M 30000 R CMD BATCH --no-save --args --\$LSB_JOBINDEX \
        --$drug --ALL genecnlc50.R log/${timestamp}_${drug}_01_02_05_\$LSB_JOBINDEX.Rout

    bsub -P $project -q ${queue} -sp ${priority} -J ${timestamp}_${drug}_01_02_06[1] \
        -oo "log/${timestamp}_${drug}_01_02_06.o" \
        -R "rusage[mem=30000]" -M 30000 -w "ended(${timestamp}_${drug}_01_02_05)" \
        $execbin/checkconcat.sh $wd/$drug/all_cn_lc50 $wd/$drug/all_cn_lc50.tsv

    bsub -P $project -q ${queue} -sp ${priority} -J ${timestamp}_${drug}_01_02_07[1] \
        -app R-2.14.0 -oo "log/${timestamp}_${drug}_01_02_07.o" \
        -w "ended(${timestamp}_${drug}_01_02_02)&&ended(${timestamp}_${drug}_01_02_04)&&ended(${timestamp}_${drug}_01_02_06)" \
        -R "rusage[mem=30000]" -M 30000 R CMD BATCH --no-save --args --$drug \
        cnlc50merge.R log/${timestamp}_${drug}_01_02_07.Rout

    bsub -P $project -q ${queue} -sp ${priority} -J ${timestamp}_${drug}_01_02_08[1-200] \
        -app R-2.14.0 -oo "log/${timestamp}_${drug}_01_02_08_%I.o" \
        -w "ended(${timestamp}_${drug}_01_02_07)" -R "rusage[mem=30000]" -M 30000 \
        R CMD BATCH --no-save --args --\$LSB_JOBINDEX --$drug \
        cnlc50meta.R log/${timestamp}_${drug}_01_02_08_\$LSB_JOBINDEX.Rout

    bsub -P $project -q ${queue} -sp ${priority} -J ${timestamp}_${drug}_01_02_09[1] \
        -oo "log/${timestamp}_${drug}_01_02_09.o" \
        -R "rusage[mem=30000]" -M 30000 -w "ended(${timestamp}_${drug}_01_02_08)" \
        $execbin/checkconcat.sh $wd/$drug/cn_lc50_meta $wd/$drug/cn_lc50.tsv

    bsub -P $project -q ${queue} -sp ${priority} -J ${timestamp}_${drug}_01_02_10[1] \
        -oo "log/${timestamp}_${drug}_01_02_10.o" -app R-2.14.0 \
        -R "rusage[mem=30000]" -M 30000  \
        -w "ended(${timestamp}_${drug}_01_02_09)" \
        R CMD BATCH --no-save --args ${drug}_cncutoff.R \
        log/${timestamp}_${drug}_01_02_10.Rout	
fi
#################################################################################

#TODO(spaugh):Reorder modules so shorter modules run first

#################################################################################
if [[ "${runge}" == "TRUE" ]]; then

    cp ${execbin}/gecutoff.R ${wd}/${drug}/${drug}_gecutoff.R

    bsub -P $project -q ${queue} -sp ${priority} -J ${timestamp}_${drug}_01_03_01[1] \
        -app R-2.14.0 -oo "log/${timestamp}_${drug}_01_03_01.o" \
        -R "rusage[mem=30000]" -M 30000 R CMD BATCH --no-save --args \
        gelc50.R log/${timestamp}_${drug}_01_03_01.Rout

    bsub -P $project -q ${queue} -sp ${priority} -J ${timestamp}_${drug}_01_03_02[1] \
        -app R-2.14.0 -oo "log/${timestamp}_${drug}_01_03_02.o" \
        -R "rusage[mem=2000]" -M 2000 -w "ended(${timestamp}_${drug}_01_03_01)" \
        R CMD BATCH --no-save --args ${drug}_gecutoff.R \
        log/${timestamp}_${drug}_01_03_02.Rout	

fi
#################################################################################

#################################################################################
if [[ "${runmir}" == "TRUE" ]]; then

    cp ${execbin}/mircutoff.R ${wd}/${drug}/${drug}_mircutoff.R

    bsub -P $project -q ${queue} -sp ${priority} -J ${timestamp}_${drug}_01_04_01[1] \
        -app R-2.14.0 -oo "log/${timestamp}_${drug}_01_04_01.o" \
        -R "rusage[mem=30000]" -M 30000 R CMD BATCH --no-save --args \
        mirlc50.R log/${timestamp}_${drug}_01_04_01.Rout

    bsub -P $project -q ${queue} -sp ${priority} -J ${timestamp}_${drug}_01_04_02[1] \
        -app R-2.14.0 -oo "log/${timestamp}_${drug}_01_04_02.o" \
        -R "rusage[mem=2000]" -M 2000 -w "ended(${timestamp}_${drug}_01_04_01)" \
        R CMD BATCH --no-save --args ${drug}_mircutoff.R \
        log/${drug}_mircutoff\$LSB_JOBINDEX.Rout

fi
#################################################################################



#################################################################################
if [[ "${runmeth}" == "TRUE" ]]; then

    cp ${execbin}/methcutoff.R ${wd}/${drug}/${drug}_methcutoff.R

    bsub -P $project -q ${queue} -sp ${priority} -J ${timestamp}_${drug}_01_05_01[1] \
        -app R-2.14.0 -oo "log/${timestamp}_${drug}_01_05_01.o" \
        -R "rusage[mem=30000]" -M 30000 R CMD BATCH --no-save --args \
        methlc50.R log/${timestamp}_${drug}_01_05_01.Rout

    bsub -P $project -q ${queue} -sp ${priority} -J ${timestamp}_${drug}_01_05_02[1] \
        -app R-2.14.0 -oo "log/${timestamp}_${drug}_01_05_02.o" \
        -R "rusage[mem=30000]" -M 30000 -w "ended(${timestamp}_${drug}_01_05_01)" \
        R CMD BATCH --no-save --args ${drug}_methcutoff.R \
        log/${timestamp}_${drug}_01_05_02.Rout

fi
#################################################################################




#################################################################################
if [[ "${rungemir}" == "TRUE" ]]; then
    mkdir -p totxv_ge_mir
    mkdir -p totxvi_ge_mir
    mkdir -p all_ge_mir
    mkdir -p ge_mir_meta
    
    cp ${execbin}/gemir.R $wd/$drug/gemir.R
    cp ${execbin}/gemirmeta.R $wd/$drug/gemirmeta.R
    
    bsub -P $project -q ${queue} -sp ${priority} -J ${timestamp}_${drug}_01_06_01[1] \
        -app R-2.14.0 -R "rusage[mem=30000]" -M 30000 \
        -oo "log/${timestamp}_${drug}_01_06_01.o" \
        -w "ended(${timestamp}_${drug}_01_03_02)&&ended(${timestamp}_${drug}_01_04_02)" \
        ${execbin}/submirge.sh ${project} ${timestamp} ${wd} ${drug} ${execbin} ${priority}
fi
#################################################################################


#################################################################################
if [[ "${rungemeth}" == "TRUE" ]]; then
    mkdir -p totxv_ge_meth
    mkdir -p totxvi_ge_meth
     mkdir -p all_ge_meth
    mkdir -p ge_meth_meta
    mkdir -p ge_meth_meta_sub
    
    cp ${execbin}/methgeprep.R ${wd}/${drug}/methgeprep.R
    cp ${execbin}/gemeth.R $wd/$drug/gemeth.R
    cp ${execbin}/gemethmeta.R $wd/$drug/gemethmeta.R
    cp ${execbin}/gemethmetasub.R $wd/$drug/gemethmetasub.R

    bsub -P $project -q ${queue} -sp ${priority} -J ${timestamp}_${drug}_7_1[1] \
        -app R-2.14.0 -oo "log/${timestamp}_${drug}_meth_ge_prep.o" \
        -R "rusage[mem=4000]" -w "ended(${timestamp}_${drug}_01_03_02)&&ended(${timestamp}_${drug}_01_05_01)" \
        R CMD BATCH --no-save --args methgeprep.R \
        log/methgeprep_${drug}_\$LSB_JOBINDEX.Rout

    bsub -P $project -q ${queue} -sp ${priority} -J ${timestamp}_${drug}_D_1[1-200] \
        -app R-2.14.0 -R "rusage[mem=4000]" -w "ended(${timestamp}_${drug}_7_1)" \
        -oo "log/${timestamp}_${drug}_totxv_meth_ge_%I.o" \
        R CMD BATCH --no-save --args --\${LSB_JOBINDEX} --TOTXV \
        gemeth.R log/totxv_meth_ge_\${LSB_JOBINDEX}.Rout

    bsub -P $project -q ${queue} -sp ${priority} -J ${timestamp}_${drug}_D_2[1-200] \
        -app R-2.14.0 -R "rusage[mem=4000]" -w "ended(${timestamp}_${drug}_7_1)" \
        -oo "log/${timestamp}_${drug}_totxvi_meth_ge_%I.o" \
        R CMD BATCH --no-save --args --\${LSB_JOBINDEX} --TOTXVI \
        gemeth.R log/totxvi_meth_ge_\${LSB_JOBINDEX}.Rout

  bsub -P $project -q ${queue} -sp ${priority} -J ${timestamp}_${drug}_D_3[1-200] \
        -app R-2.14.0 -R "rusage[mem=4000]" -w "ended(${timestamp}_${drug}_7_1)" \
        -oo "log/${timestamp}_${drug}_all_meth_ge_%I.o" \
        R CMD BATCH --no-save --args --\${LSB_JOBINDEX} --ALL \
        gemeth.R log/all_meth_ge_\${LSB_JOBINDEX}.Rout


    bsub -P $project -q ${queue} -sp ${priority} -J ${timestamp}_${drug}_D_4[1-200] \
        -app R-2.14.0 -R "rusage[mem=4000]" -w "ended(${timestamp}_${drug}_D_1)&&ended(${timestamp}_${drug}_D_2)&&ended(${timestamp}_${drug}_D_3)" \
        -oo "log/${timestamp}_${drug}_meth_ge_meta_%I.o" \
        R CMD BATCH --no-save --args --\${LSB_JOBINDEX} gemethmeta.R \
        log/meth_ge_meta_\${LSB_JOBINDEX}.Rout

    bsub -P $project -q ${queue} -sp ${priority} -J ${timestamp}_${drug}_D_5[1-200] \
        -app R-2.14.0 -R "rusage[mem=4000]" -w "ended(${timestamp}_${drug}_D_4)" \
        -oo "log/${timestamp}_${drug}_meth_ge_meta_sub_%I.o" \
        R CMD BATCH --no-save --args --\${LSB_JOBINDEX} gemethmetasub.R \
        log/meth_ge_meta_sub_${drug}_\${LSB_JOBINDEX}.Rout

    bsub -P $project -q ${queue} -sp ${priority} -J ${timestamp}_${drug}_D_6[1] \
        -w "ended(${timestamp}_${drug}_D_5)" -R "rusage[mem=30000]" \
        -oo "log/${timestamp}_${drug}_meth_ge_meta_sub_result_cat_%I.o" \
        $execbin/checkconcat.sh ${wd}/${drug}/ge_meth_meta_sub ${wd}/${drug}/ge_meth_meta_sub.tsv
    

fi
#################################################################################


#################################################################################
if [[ "${rungesnp}" == "TRUE" ]]; then
    cp ${execbin}/snpgeprep.R $wd/$drug/snpgeprep.R

    bsub -P $project -q ${queue} -sp ${priority} -J ${timestamp}_${drug}_4_1[1] \
        -app R-2.14.0 -o "log/${timestamp}_${drug}_snp_ge_prep_%I.o" \
        -w "ended(${timestamp}_${drug}_01_03_02)&&ended(${timestamp}_${drug}_01_01_14)" \
        -R "rusage[mem=30000]" R CMD BATCH --no-save --args snpgeprep.R \
        log/snpgeprep_${drug}_\$LSB_JOBINDEX.Rout

    mkdir -p totxv_ge_snp
    mkdir -p totxvi_ge_snp
    mkdir -p all_ge_snp
    mkdir -p ge_snp_meta
    mkdir -p ge_snp_meta_sub

    cp ${execbin}/gesnp.R $wd/$drug/gesnp.R
    cp ${execbin}/gesnpmeta.R $wd/$drug/gesnpmeta.R
    cp ${execbin}/gesnpmetasub.R $wd/$drug/gesnpmetasub.R

    bsub -P $project -q ${queue} -sp ${priority} -J ${timestamp}_${drug}_4_2[1] \
        -app R-2.14.0 -R "rusage[mem=30000]" -w "ended(${timestamp}_${drug}_4_1)" \
        ${execbin}/subsnpge.sh ${project} ${timestamp} ${wd} ${drug} ${execbin} ${queue}  ${priority}
fi
#################################################################################


#################################################################################
if [[ "${rungecn}" == "TRUE" ]]; then

    cp ${execbin}/cngeprep.R $wd/$drug/cngeprep.R

    bsub -P $project -q ${queue} -sp ${priority} -J ${timestamp}_${drug}_5_1[1] \
        -w "ended(${timestamp}_${drug}_01_03_02)&&ended(${timestamp}_${drug}_01_02_10)" \
        -app R-2.14.0 -oo "log/${timestamp}_${drug}_cn_ge_prep_%I.o" \
        -R "rusage[mem=30000]" R CMD BATCH --no-save --args \
        cngeprep.R log/cngeprep_${drug}_\$LSB_JOBINDEX.Rout

    mkdir -p totxv_ge_cn
    mkdir -p totxvi_ge_cn
    mkdir -p all_ge_cn
    mkdir -p ge_cn_meta
    mkdir -p ge_cn_meta_sub

    cp ${execbin}/genegecn.R $wd/$drug/gecn.R
    cp ${execbin}/gecnmeta.R $wd/$drug/gecnmeta.R
    cp ${execbin}/gecnmetasub.R $wd/$drug/gecnmetasub.R

    bsub -P $project -q ${queue} -sp ${priority} -J ${timestamp}_${drug}_5_2[1] \
        -app R-2.14.0 -R "rusage[mem=10000]" \
        -w "ended(${timestamp}_${drug}_5_1)" \
        ${execbin}/subcnge.sh ${project} ${timestamp} ${wd} ${drug} ${execbin} ${queue}

fi
#################################################################################

#################################################################################
if [[ "${runmirsnp}" == "TRUE" ]]; then

    cp ${execbin}/snpmirprep.R $wd/$drug/snpmirprep.R

   # bsub -P $project -q ${queue} -J ${timestamp}_${drug}_8_1[1] -app R-2.14.0 -o "log/${timestamp}_${drug}_snp_mir_prep_%I.o" -R "rusage[mem=30000]" R CMD BATCH --no-save --args snpmirprep.R log/snpmirprep_${drug}_\$LSB_JOBINDEX.Rout

    mkdir -p totxv_snp_mir
    mkdir -p totxvi_snp_mir
    mkdir -p snp_mir_meta
    mkdir -p snp_mir_meta_sub

    cp ${execbin}/snpmir.R $wd/$drug/snpmir.R
    cp ${execbin}/snpmirmeta.R $wd/$drug/snpmirmeta.R
    cp ${execbin}/snpmirmetasub.R $wd/$drug/snpmirmetasub.R

  # bsub -P $project -q ${queue} -J ${timestamp}_${drug}_8_2[1] -app R-2.14.0 -R "rusage[mem=30000]" -w "ended(${timestamp}_${drug}_8_1)" ${execbin}/subsnpmir.sh ${project} ${timestamp} ${wd} ${drug} ${execbin} ${queue}

fi
#################################################################################



#################################################################################
if [[ "${runfigure}" == "TRUE" ]]; then
    
    bsub -P $project -q ${queue} -sp ${priority} -J ${timestamp}_${drug}_7_1[1] \
        -app R-2.14.0 -R "rusage[mem=1000]" -M 1000 \
        -oo "log/${timestamp}_${drug}_01_07_01.o" \
        -w "done(${timestamp}_${drug}_01_01_14)&&done(${timestamp}_${drug}_01_02_10)&&done(${timestamp}_${drug}_01_03_02)&&done(${timestamp}_${drug}_01_04_02)&&done(${timestamp}_${drug}_01_05_02)&&done(${timestamp}_${drug}_01_06_01)" \
        ${execbin}/subfigure.sh ${project} ${timestamp} ${wd} ${drug} ${execbin} ${priority}
fi
#################################################################################

#################################################################################
if [[ "${runonlyfigure}" == "TRUE" ]]; then

    bsub -P $project -q ${queue} -sp ${priority} -J ${timestamp}_${drug}_7_1[1] \
        -app R-2.14.0 -R "rusage[mem=1000]" -M 1000 \
        -oo "log/${timestamp}_${drug}_01_07_01.o" \
        ${execbin}/subfigure.sh ${project} ${timestamp} ${wd} ${drug} ${execbin} ${priority}

fi
#################################################################################
#################################################################################
if [[ "${runbabyfigure}" == "TRUE" ]]; then
    
    bsub -P $project -q ${queue} -sp ${priority} -J ${timestamp}_${drug}_8_1[1] \
        -app R-2.14.0 -R "rusage[mem=1000]" -M 1000 \
        -oo "log/${timestamp}_${drug}_01_08_01.o" \
        -w "ended(${timestamp}_${drug}_01_01_14)&&ended(${timestamp}_${drug}_01_02_10)&&ended(${timestamp}_${drug}_01_03_02)&&ended(${timestamp}_${drug}_01_04_02)&&ended(${timestamp}_${drug}_01_05_02)&&ended(${timestamp}_${drug}_01_06_01)&&ended(${timestamp}_${drug}_01_07_01)" \
        ${execbin}/subbabyfigure.sh ${project} ${timestamp} ${wd} ${drug} ${execbin} ${priority}
fi
#################################################################################

#################################################################################
if [[ "${runonlybabyfigure}" == "TRUE" ]]; then

    bsub -P $project -q ${queue} -sp ${priority} -J ${timestamp}_${drug}_8_1[1] \
        -app R-2.14.0 -R "rusage[mem=1000]" -M 1000 \
        -oo "log/${timestamp}_${drug}_01_08_01.o" \
        ${execbin}/subbabyfigure.sh ${project} ${timestamp} ${wd} ${drug} ${execbin} ${priority}

fi
################################################################################
#################################################################################
if [[ "${runpackageresult}" == "TRUE" ]]; then

echo ""
#bsub -P $project -q ${queue} -J ${timestamp}_${drug}_resultprep_1 -oo "log/${timestamp}_${drug}_resultprep_1.o" -R "rusage[mem=2000]" ${execbin}/result_prep.sh ${wd} ${drug} ${timestamp}

fi
#################################################################################

done
