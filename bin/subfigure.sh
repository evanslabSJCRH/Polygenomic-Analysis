#!/bin/bash

project=${1}
timestamp=${2}
wd=${3}
drug=${4}
execbin=${5}
priority=${6}

queue="standard"
lmqueue="standard"
#lmhost="uv02"
#TODO(spaugh): Source config.example

#################################################################################
while  [  ! -f ge_mir_meta.tsv ]
do
    sleep 100
done

while  [  ! -f ge_meth_meta_sub.tsv ]
do
    sleep 100
done

#cp ${execbin}/tablecreate.Rnw $wd/$drug/tablecreate.Rnw

cd ${wd}/${drug}

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


create_sweave_dispatch tablecreate

bsub -P $project -sp ${priority} -q ${queue} -J ${timestamp}_${drug}_table_01[1] \
    -app R-2.14.0 -oo "log/${timestamp}_${drug}_table_01_%I.o" \
    -R "rusage[mem=4000]" -M 4000 R CMD BATCH --no-save --args tablecreate.R \
    log/table_\$LSB_JOBINDEX.Rout

bsub -P $project -sp ${priority} -q ${queue} -J ${timestamp}_${drug}_table_pdf_1[1] \
    -app R-2.14.0 -oo "log/${timestamp}_${drug}_table_pdf_%I.o" \
    -R "rusage[mem=1000]" -M 1000 -w "ended(${timestamp}_${drug}_table_01)" \
    pdflatex tablecreate.tex	

bsub -P $project -sp ${priority} -q ${queue} -J ${timestamp}_${drug}_table_pdf_2[1] \
    -app R-2.14.0 -oo "log/${timestamp}_${drug}_table_pdf_%I.o" \
    -R "rusage[mem=1000]" -M 1000 -w "ended(${timestamp}_${drug}_table_pdf_1)" \
    pdflatex tablecreate.tex	

bsub -P $project -sp ${priority} -q ${queue} -J ${timestamp}_${drug}_table_pdf_3[1] \
    -app R-2.14.0 -oo "log/${timestamp}_${drug}_table_pdf_%I.o" \
    -R "rusage[mem=1000]" -M 1000 -w "ended(${timestamp}_${drug}_table_pdf_2)" \
    pdflatex tablecreate.tex	
#################################################################################

cp ${execbin}/figures.Rnw ${wd}/${drug}/figures_${drug}.Rnw

create_sweave_dispatch figures
    
    bsub -P $project -sp ${priority} -q ${lmqueue} -J ${timestamp}_${drug}_figures[1] \
        -oo "figures.o" -app R-2.14.0 -R "rusage[mem=100000]" -M 100000 \
        R CMD BATCH --no-save --args figures.R \
        log/${drug}_figures_\$LSB_JOBINDEX.Rout


#################################################################################
