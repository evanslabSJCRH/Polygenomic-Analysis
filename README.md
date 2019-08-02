
This repository contains code to analyze LC50 data against a variety of different genomic features.

## Design

The design of the analysis breaks each genomic feature or connection between genomic features into distinct modules.  These modules are referenced in the job name and used as completion dependencies throughout the workflow.  The workflow is separated into three phases.  In the first phase, LC50 is analyzed against the genomic features and features are selected.  In the second phase, these features are tested for connectivity. In the third phase, figures are drawn and analyses are summarized.

## Style

When possible the Google style guide for code is used.  There are exceptions but their style guides provide a good baseline for which to strive.

Links to pertinent guides:

[Bash](http://google-styleguide.googlecode.com/svn/trunk/shell.xml)  
Notable exceptions:  
Use 4 space indents rather than 2 space indents.  
Prefered file extension is .sh rather than no file extension.  

[R](http://google-styleguide.googlecode.com/svn/trunk/Rguide.xml)

## Running the workflow

### Initial setup

To execute the workflow you first need to clone this repository.



### Execution of workflow

Run the workflow using the following command:

    ./run02.sh

## Workflow

### 01. Phase 01

#### Setup

The begining of the workflow setups the output directories, some bash variables and copies or symlinks files into the output directories.

Some notable variables setup and read from the config file are:

    ${wd}:   This is the working directory where results are output
    ${drug}: This is the current drug being run 

#### 01.01. SNP vs LC50

##### 01.01.01 TOTXV SNP vs LC50

*Dependencies: None*

This section uses the snplc50.R script.

This will create 200 files with the results of the GWAS for SNPs versus LC50.  These files will be placed into the totxv_snp_lc50 folder with filenames output.001 - output.200.  Each file has a header row and is tab delimited.

Example output:

    ProbeSetID     fisher.p.b         fisher.01.b  fisher.03.b  fisher.11.b  fisher.13.b  fisher.21.b  fisher.23.b  snp.lm.p.b         snp.lm.stat.b
    SNP_A-8593276  0.5684297996627    54           17           4            0            0            0            0.418536020081667  -0.812500320587933
    SNP_A-8593277  0.568494131711743  1            1            10           3            47           13           0.539813706357981  -0.615321458131022
    SNP_A-8593278  0.32068163538213   3            0            20           3            34           13           0.0828162011995449  1.75367728333007
    SNP_A-8593279  0.559018748198693  2            1            5            2            50           14           0.389910553985738  -0.863774395168782

The columns represent the following results:

    ProbeSetID:    This is the Affymetrix ProbeSetID
    fisher.p.b:    This is the Fisher test p-value
    fisher.01.b:   This is the number of patients that are homozygous   AA (0) and in the drug sensitive group (1)
    fisher.03.b:   This is the number of patients that are homozygous   AA (0) and in the drug resistant group (3)
    fisher.11.b:   This is the number of patients that are heterozygous AB (1) and in the drug sensitive group (1)
    fisher.13.b:   This is the number of patients that are heterozygous AB (1) and in the drug resistant group (3)
    fisher.21.b:   This is the number of patients that are homozygous   BB (2) and in the drug sensitive group (1)
    fisher.23.b:   This is the number of patients that are homozygous   BB (2) and in the drug resistant group (3)
    snp.lm.p.b:    This is the linear model p-value
    snp.lm.stat.b: This is the linear model test statistic

##### 01.01.02 TOTXV SNP vs LC50 result concatentation

*Dependencies: Section 01.01.01 completion*

This section uses the general purpose concat.sh script which takes all files in a folder, and combines them into a single file removing extranous headers from all but the first file.  This uses the folder **${wd}/${drug}/totxv_snp_lc50** as input and outputs a single file **${wd}/${drug}/totxv_snp_lc50.tsv**.

##### 01.01.03 TOTXVI SNP vs LC50

*Dependencies: None*

This section is similar to section 01.01.01 and uses the same snplc50.R script.  Output files are placed into the **${wd}/${drug}/totxvi_snp_lc50** folder similar to section 01.01.01.

##### 01.01.04 TOTXVI SNP vs LC50 result concatentation

*Dependencies: Section 01.01.03 completion*

This section is similar to section 01.01.02 and uses the same concat.sh script.  This uses the folder **${wd}/${drug}/totxvi_snp_lc50** as input and outputs a single file **${wd}/${drug}/totxvi_snp_lc50.tsv**.

##### 01.01.05 ALL patients SNP vs LC50

*Dependencies: None*

This section is similar to section 01.01.01 and uses the same snplc50.R script.  Output files are placed into the **${wd}/${drug}/all_snp_lc50** folder similar to section 01.01.01.  Patients from both TOTXV and TOTXVI are used in one analysis using protocol as a covariate.

##### 01.01.06 ALL patients SNP vs LC50 result concatentation

*Dependencies: Section 01.01.05 completion*

This section is similar to section 01.01.02 and uses the same concat.sh script.  This uses the folder **${wd}/${drug}/all_snp_lc50** as input and outputs a single file **${wd}/${drug}/all_snp_lc50.tsv**.

##### 01.01.07 Imputted TOTXV SNP vs LC50

*Dependencies: None*

This section is similar to section 01.01.01 and uses the same snplc50.R script.  Output files are placed into the **${wd}/${drug}/imputedxv_snp_lc50** folder similar to section 01.01.01.

##### 01.01.08 Imputted TOTXV SNP vs LC50 result concatentation

##### 01.01.09 Imputted TOTXVI SNP vs LC50

*Dependencies: None*

This section is similar to section 01.01.01 and uses the same snplc50.R script.  Output files are placed into the **${wd}/${drug}/imputedxv_snp_lc50** folder similar to section 01.01.01.

##### 01.01.10 Imputted TOTXVI SNP vs LC50 result concatentation

##### 01.01.11 SNP versus LC50 result merge

*Dependencies: Sections 01.01.02, 01.01.04, 01.01.06 completion*

##### 01.01.12 TOTXV, TOTXVI meta analysis
##### 01.01.13 TOTXV, TOTXVI meta analysis result concatentation
##### 01.01.14 SNP cutoff


#### 01.02. CNV vs LC50


##### 01.02.01 TOTXV CNV vs LC50

*Dependencies: None*

This section uses the cnlc50.R script.

This will create 200 files with the results of the GWAS for CNs versus LC50.  These files will be placed into the totxv_cn_lc50 folder with filenames output.001 - output.200.  Each file has a header row and is tab delimited.

    Probe.Set.ID  t.p.b            t.stat.b          w.p.b             w.stat.b          cn.lm.p.b         cn.lm.stat.b
    CN_473963     0.697868280636   0.390108930714258 0.836626263578179 0.206814846585984 0.687634707600767 0.403298769918359
    CN_473964     0.77729540526443 0.28421742899216  0.853492484193314 0.196598861613929 0.711721604222232 0.37065244117726
    CN_473965     0.77729540526443 0.28421742899216  0.853492484193314 0.196598861613929 0.711721604222232 0.37065244117726
    CN_473981     0.77729540526443 0.28421742899216  0.853492484193314 0.196598861613929 0.711721604222232 0.37065244117726
    CN_473982     0.77729540526443 0.28421742899216  0.853492484193314 0.196598861613929 0.711721604222232 0.37065244117726



##### 01.02.02 TOTXV CNV vs LC50 result concatentation


##### 01.02.03 TOTXVI CNV vs LC50
##### 01.02.04 TOTXVI CNV vs LC50 result concatentation
##### 01.02.05 ALL CNV vs LC50
##### 01.02.06 ALL CNV vs LC50 result concatentation
##### 01.02.07 TOTXV, TOTXVI, ALL result merge
##### 01.02.08 TOTXV, TOTXVI meta analysis
##### 01.02.09 TOTXV, TOTXVI meta analysis result concatentation
##### 01.02.10 CNV cutoff

#### 01.03. Gene expression vs LC50

##### 01.03.01 GE vs LC50

*Dependencies: None*

This section uses the gelc50.R and gelc50.Rnw scripts.  gelc50.R is created when executing run01.sh.  gelc50.Rnw has the analysis code.  These Rnw files are used in conjection with the R files so that result caching may be used.

This will create output files with the results of the GWAS for gene expression versus LC50.  These files will be placed into the root of the drug output folder (e.g. "/home/rautry/workflow_results_final/PRED/" ).  An output file is created for each cohort analyzed and a combined file.  

Output files:

    dutch_ge_lc50.tsv
    totxv_ge_lc50.tsv
    totxvi_ge_lc50.tsv
    ge_lc50.tsv 
    

Output file columns (dutch_ge_lc50.tsv):

    Probe.Set.ID
    dutch.t.p.bt
    dutch.t.stat.bt
    dutch.w.p.bt
    dutch.w.stat.bt
    dutch.lm.p.bt
    dutch.lm.stat.bt
    dutch.t.p.b
    dutch.t.stat.b
    dutch.w.p.b
    dutch.w.stat.b
    dutch.lm.p.b
    dutch.lm.stat.b
    dutch.t.p.t
    dutch.t.stat.t
    dutch.w.p.t
    dutch.w.stat.t
    dutch.lm.p.t
    dutch.lm.stat.t

Output file columns (totxv_ge_lc50.tsv):

    Probe.Set.ID
    totxv.t.p.bt
    totxv.t.stat.bt
    totxv.w.p.bt
    totxv.w.stat.bt
    totxv.lm.p.bt
    totxv.lm.stat.bt
    totxv.t.p.b
    totxv.t.stat.b
    totxv.w.p.b
    totxv.w.stat.b
    totxv.lm.p.b
    totxv.lm.stat.b
    totxv.t.p.t
    totxv.t.stat.t
    totxv.w.p.t
    totxv.w.stat.t
    totxv.lm.p.t
    totxv.lm.stat.t

Output file columns (totxvi_ge_lc50.tsv):

    Probe.Set.ID
    totxvi.t.p.bt
    totxvi.t.stat.bt
    totxvi.w.p.bt
    totxvi.w.stat.bt
    totxvi.lm.p.bt
    totxvi.lm.stat.bt
    totxvi.t.p.b
    totxvi.t.stat.b
    totxvi.w.p.b
    totxvi.w.stat.b
    totxvi.lm.p.b
    totxvi.lm.stat.b
    totxvi.t.p.t
    totxvi.t.stat.t
    totxvi.w.p.t
    totxvi.w.stat.t
    totxvi.lm.p.t
    totxvi.lm.stat.t

Output file columns (ge_lc50.tsv):

    Probe.Set.ID
    totxv.t.p.bt
    totxv.t.stat.bt
    totxv.w.p.bt
    totxv.w.stat.bt
    totxv.lm.p.bt
    totxv.lm.stat.bt
    totxv.t.p.b
    totxv.t.stat.b
    totxv.w.p.b
    totxv.w.stat.b
    totxv.lm.p.b
    totxv.lm.stat.b
    totxv.t.p.t
    totxv.t.stat.t
    totxv.w.p.t
    totxv.w.stat.t
    totxv.lm.p.t
    totxv.lm.stat.t
    totxvi.t.p.bt
    totxvi.t.stat.bt
    totxvi.w.p.bt
    totxvi.w.stat.bt
    totxvi.lm.p.bt
    totxvi.lm.stat.bt
    totxvi.t.p.b
    totxvi.t.stat.b
    totxvi.w.p.b
    totxvi.w.stat.b
    totxvi.lm.p.b
    totxvi.lm.stat.b
    totxvi.t.p.t
    totxvi.t.stat.t
    totxvi.w.p.t
    totxvi.w.stat.t
    totxvi.lm.p.t
    totxvi.lm.stat.t
    all.lm.p.bt
    all.lm.stat.bt
    sj.lm.p.bt
    sj.lm.stat.bt
    all.lm.p.b
    all.lm.stat.b
    sj.lm.p.b
    sj.lm.stat.b
    all.lm.p.t
    all.lm.stat.t
    sj.lm.p.t
    sj.lm.stat.t
    meta.w.stat.bt
    meta.w.p.bt
    meta.w.stat.b
    meta.w.p.b
    meta.w.stat.t
    meta.w.p.t
    meta.t.stat.bt
    meta.t.p.bt
    meta.t.stat.b
    meta.t.p.b
    meta.t.stat.t
    meta.t.p.t
    meta.lm.stat.bt
    meta.lm.p.bt
    meta.lm.stat.b
    meta.lm.p.b
    meta.lm.stat.t
    meta.lm.p.t
    stat.bt
    p.bt
    stat.b
    p.b
    stat.t
    p.t
    

##### 01.03.02 GE cutoff

*Dependencies: Section 3.02 completion*

#### 01.04. MicroRNA vs LC50

##### 01.04.01 MicroRNA vs LC50

*Dependencies: None*

Output file columns (totxv_mir_lc50.tsv):

    Name
    totxv.mir.lm.p.bt
    totxv.mir.lm.stat.bt
    totxv.mir.t.p.bt
    totxv.mir.t.stat.bt
    totxv.mir.w.p.bt
    totxv.mir.w.stat.bt
    totxv.max.p.bt
    totxv.mir.lm.p.b
    totxv.mir.lm.stat.b
    totxv.mir.t.p.b
    totxv.mir.t.stat.b
    totxv.mir.w.p.b
    totxv.mir.w.stat.b
    totxv.max.p.b
    totxv.mir.lm.p.t
    totxv.mir.lm.stat.t
    totxv.mir.t.p.t
    totxv.mir.t.stat.t
    totxv.mir.w.p.t
    totxv.mir.w.stat.t
    totxv.max.p.t
    

Output file columns (totxvi_mir_lc50.tsv):

    Name
    totxvi.mir.lm.p.bt
    totxvi.mir.lm.stat.bt
    totxvi.mir.t.p.bt
    totxvi.mir.t.stat.bt
    totxvi.mir.w.p.bt
    totxvi.mir.w.stat.bt
    totxvi.max.p.bt
    totxvi.mir.lm.p.b
    totxvi.mir.lm.stat.b
    totxvi.mir.t.p.b
    totxvi.mir.t.stat.b
    totxvi.mir.w.p.b
    totxvi.mir.w.stat.b
    totxvi.max.p.b
    totxvi.mir.lm.p.t
    totxvi.mir.lm.stat.t
    totxvi.mir.t.p.t
    totxvi.mir.t.stat.t
    totxvi.mir.w.p.t
    totxvi.mir.w.stat.t
    totxvi.max.p.t
        
Output file columns (mir_lc50.tsv):

    Name
    totxv.mir.lm.p.bt
    totxv.mir.lm.stat.bt
    totxv.mir.t.p.bt
    totxv.mir.t.stat.bt
    totxv.mir.w.p.bt
    totxv.mir.w.stat.bt
    totxv.max.p.bt
    totxv.mir.lm.p.b
    totxv.mir.lm.stat.b
    totxv.mir.t.p.b
    totxv.mir.t.stat.b
    totxv.mir.w.p.b
    totxv.mir.w.stat.b
    totxv.max.p.b
    totxv.mir.lm.p.t
    totxv.mir.lm.stat.t
    totxv.mir.t.p.t
    totxv.mir.t.stat.t
    totxv.mir.w.p.t
    totxv.mir.w.stat.t
    totxv.max.p.t
    totxvi.mir.lm.p.bt
    totxvi.mir.lm.stat.bt
    totxvi.mir.t.p.bt
    totxvi.mir.t.stat.bt
    totxvi.mir.w.p.bt
    totxvi.mir.w.stat.bt
    totxvi.max.p.bt
    totxvi.mir.lm.p.b
    totxvi.mir.lm.stat.b
    totxvi.mir.t.p.b
    totxvi.mir.t.stat.b
    totxvi.mir.w.p.b
    totxvi.mir.w.stat.b
    totxvi.max.p.b
    totxvi.mir.lm.p.t
    totxvi.mir.lm.stat.t
    totxvi.mir.t.p.t
    totxvi.mir.t.stat.t
    totxvi.mir.w.p.t
    totxvi.mir.w.stat.t
    totxvi.max.p.t
    all.mir.lm.p.bt
    all.mir.lm.stat.bt
    all.mir.lm.p.b
    all.mir.lm.stat.b
    all.mir.lm.p.t
    all.mir.lm.stat.t
    meta.stat.bt
    meta.p.bt
    meta.stat.b
    meta.p.b
    meta.stat.t
    meta.p.t
    p.bt
    stat.bt
    p.b
    stat.b
    p.t
    stat.t
    
##### 01.04.02 MicroRNA cutoff

#### 01.05. Methylation vs LC50

##### 01.05.01 Methylation vs LC50
*Dependencies: None*

For the catagorical tests and TOTXV DEX 6MP LASP VCR T-only there not enough observations so we do not do these tests. Use of catagorical tests is deprecated, but results are still output for some portions.   Reporting and testing will be eliminated for catagorical tests.

##### 01.05.02 Methylation cutoff
*Dependencies: Section 5.01 completion*


### Phase 2

#### 4. SNP vs GE
##### 4.1 SNP vs GE prep
##### 4.2 Submission of dispatch script

#### 5. CN vs GE
##### 5.1 CN vs GE prep
##### 5.2 Submission of dispatch script

#### 6. MIR vs GE
##### 6.1 MIR vs GE prep
##### 6.2 Submission of dispatch script 


#### 7. METH vs GE
##### 7.1 METH vs GE prep

#### D. METH vs GE
##### D.1 TOTXV METH vs GE prep
##### D.2 TOTXVI METH vs GE prep
##### D.3 METH vs GE meta
##### D.4 METH vs GE meta subset
##### D.5 METH vs GE meta subset concatentation

#### 8. SNP vs MIR
#### 8.1 SNP vs MIR prep
#### 8.2 Submission of dispatch script 


### Phase 3


## Monitoring job progress

At the beginning of each phase dispatch a timestamp is echoed which is used in the construction of the job names.

To monitor job status "bjobs -A -J timestamp*" can be used to get a high level overview of indexed job status.

Example dispatch output from Phase 1:

    [rautry@erus02 drworkflow]$ ./run01.sh 
    /home/rautry/workflow_results_final
    20131114200328874045195
    PRED
    Job <6438573> is submitted to queue <priority>.
    Job <6438574> is submitted to queue <priority>.
    Job <6438575> is submitted to queue <priority>.
    Job <6438576> is submitted to queue <priority>.
    Job <6438577> is submitted to queue <priority>.
    Job <6438578> is submitted to queue <priority>.
    Job <6438579> is submitted to queue <priority>.
    Job <6438580> is submitted to queue <priority>.
    Job <6438581> is submitted to queue <priority>.
    Job <6438582> is submitted to queue <priority>.
    Job <6438583> is submitted to queue <priority>.
    Job <6438584> is submitted to queue <priority>.
    Job <6438585> is submitted to queue <priority>.
    Job <6438586> is submitted to queue <priority>.
    Job <6438587> is submitted to queue <priority>.
    Job <6438588> is submitted to queue <priority>.
    Job <6438589> is submitted to queue <priority>.
    Job <6438590> is submitted to queue <priority>.
    Job <6438591> is submitted to queue <priority>.
    Job <6438592> is submitted to queue <priority>.
    Job <6438593> is submitted to queue <priority>.
    Job <6438594> is submitted to queue <priority>.
    Job <6438595> is submitted to queue <priority>.
    Job <6438596> is submitted to queue <priority>.
    Job <6438597> is submitted to queue <uv>.
    
Example monitoring of Phase 1 dispatch:

    [rautry@erus02 drworkflow]$ bjobs -A -J 20131114200328874045195*
    JOBID       ARRAY_SPEC  OWNER   NJOBS PEND DONE  RUN EXIT SSUSP USUSP PSUSP
    6438573     20131114   rautry     200  104   95    1    0     0     0     0
    6438574     20131114   rautry       1    1    0    0    0     0     0     0
    6438575     20131114   rautry     200  200    0    0    0     0     0     0
    6438576     20131114   rautry       1    1    0    0    0     0     0     0
    6438577     20131114   rautry     200  200    0    0    0     0     0     0
    6438578     20131114   rautry       1    1    0    0    0     0     0     0
    6438579     20131114   rautry       1    1    0    0    0     0     0     0
    6438580     20131114   rautry     200  200    0    0    0     0     0     0
    6438581     20131114   rautry       1    1    0    0    0     0     0     0
    6438582     20131114   rautry     200  200    0    0    0     0     0     0
    6438583     20131114   rautry       1    1    0    0    0     0     0     0
    6438584     20131114   rautry     200  200    0    0    0     0     0     0
    6438585     20131114   rautry       1    1    0    0    0     0     0     0
    6438586     20131114   rautry     200  200    0    0    0     0     0     0
    6438587     20131114   rautry       1    1    0    0    0     0     0     0
    6438588     20131114   rautry       1    1    0    0    0     0     0     0
    6438589     20131114   rautry     200  200    0    0    0     0     0     0
    6438590     20131114   rautry       1    1    0    0    0     0     0     0
    6438591     20131114   rautry       1    1    0    0    0     0     0     0
    6438592     20131114   rautry       1    1    0    0    0     0     0     0
    6438593     20131114   rautry       1    1    0    0    0     0     0     0
    6438594     20131114   rautry       1    1    0    0    0     0     0     0
    6438595     20131114   rautry       1    1    0    0    0     0     0     0
    6438596     20131114   rautry       1    1    0    0    0     0     0     0
    6438597     20131114   rautry       1    1    0    0    0     0     0     0

## Working with TSV files

TSV files are Tab Separated Value files.  In short, each field is separated by a tab.  TSV files are essentially text files and can be viewed by a wide variety of programs.

To read and write in the "preferred" TSV format using R see the example lines of code below.

    write.table (output.result, file="output_result.tsv", row.names=FALSE, quote=FALSE, sep="\t")
    output.result <- read.delim ("output_result.tsv", as.is=TRUE, stringsAsFactors=FALSE)

## Output file naming convention

Output files are useful tools for troubleshooting.  Output files generated by bsub should end with a ".o" and be named similar to the following example:

    "log/${timestamp}_${drug}_01_01_01_%I.o"
    "log/${timestamp}_${drug}_01_01_02.o"

Where the numbers indicates the Phase, Section and Subsection number respectively, with leading zeros for numbers less than 10.  We should not to have more than two digits for the phase, section or subsection number. 
For jobs with a single index, *omit* the trailing index identifier.  For jobs with multiple indexes, include the index identifier with a "_%I".