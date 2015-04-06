import sys
import subprocess
import numpy as np
import csv
import os

def main():
    
    print "Hello World\n"

    #FLAGS
    flags = {};
    setDefaults(flags);
    
    #Input Args
    (flags, projName, sampleName, f1, f2, logfileName, sampleNum) = parseInput(sys.argv, flags)
    
    #Directories
    (baseBigDir, samplesDir,bamDir, tdfDir, cuffDir, logDir, metricsDir, jobDir, runAllDir, RSEM_AllOutDir, Cuff_AllOutDir) = createDirs(flags,projName)
   
    #Overall log
    if sampleNum==1:
        logfile = open(baseBigDir + '/' + logfileName, "w+");
    else: 	
        logfile = open(baseBigDir + '/' + logfileName, "a+");
  	
    logWrite(logfile,"Running pipeline for sample " + sampleName)
    
    #Sample Log
    SampleLog = open(logDir + '/' + 'log_' + sampleName + '.txt',"w+")
    
    #Create job file
    jobfile = open(jobDir + '/' + 'job_' + sampleName + '.sh', "w+")
    JobFileHeader(jobfile)

    logWrite(SampleLog,"Starting pipeline for sample " + sampleName);
    
    #FASTq extraction

    #STARTING PIPELINE
    #If input is fastq.gz, then unzip the file
    if '.gz' in f1:
        gzipFlag = 1
	logWrite(jobfile,"\n#Unzipping fastq files for sample " + sampleName + "\n");
        if flags['single-end']==1:
	    logWrite(jobfile, "gunzip " + f1)
	    f1 = f1[:-3]
	else:
	    logWrite(jobfile, "gunzip " + f1)
	    logWrite(jobfile, "gunzip " + f2)
	    f1 = f1[:-3]; f2 = f2[:-3]
    else:
        gzipFlag = 0

    #Confirm existence of fastq files

    #SMART trim

    if flags['doTophat'] == '1':
        logWrite(SampleLog,"Running Tophat for sample " + sampleName);
        (finalBam, finalBamLink, oldBam) = runTophat(jobfile, flags, sampleName, samplesDir, bamDir, fq1=f1, fq2=f2)
        
	#save tophat results
    if flags['doPicard'] == '1':
        logWrite(SampleLog, "Running Picard tools on Tophat output for sample " + sampleName);
        runPicard(jobfile, flags, sampleName, metricsDir, bamName=finalBamLink)	

    if flags['doBam2Tdf'] == '1':
        logWrite(SampleLog,"Converting BAM to TDF");
	Bam2Tdf(jobfile, flags, sampleName, tdfDir, bamName=finalBamLink)

    if flags['doRSEM'] == '1':
        logWrite(SampleLog,"Running RSEM for sample " + sampleName);
        runRSEM(jobfile, flags, sampleName, samplesDir,metricsDir, RSEM_AllOutDir, fq1=f1, fq2=f2)

    if flags['doCufflinks'] == '1':
        logWrite(SampleLog,"Running Cufflinks for sample " + sampleName);
	runCufflinks(jobfile, flags, sampleName, samplesDir, Cuff_AllOutDir, bamName=finalBamLink)

    if flags['doRNASeqQC'] == '1':
        logWrite(SampleLog,"Running RNASeqQC for sample " + sampleName);
	runRNASeqQC(jobfile, flags, sampleName, metricsDir, bamName=finalBamLink)
    
    if flags['doCountrRNA'] == '1':
        logWrite(SampleLog, "Checking for rRNA in sample " + sampleName);
	runCountrRNA(jobfile, flags, sampleName, metricsDir, fq1=f1, fq2=f2)
    if flags['doRSeQC'] == '1':
        logWrite(SampleLog, "Running RSeQC in sample " + sampleName);
	runRSeQC(jobfile, flags, sampleName, metricsDir, bamName=finalBamLink)

    #return
    #if flags['doBowtie'] == 'adadsad':
        #runBowtie()
	#if not os.path.exists(rsemOutDir): os.makedirs(rsemOutDir);
	#Check if all the fastq files in the list exist (TO DO)
	#logWrite(runfile, "\n Running Bowtie for RSEM \n");
    
    if gzipFlag==1:
        logWrite(jobfile, "\n\n#Re-zipping fastq files for sample " + sampleName + "\n");
        if flags["single-end"]==1:
	    logWrite(jobfile, "gzip " + f1)
        else:
	    logWrite(jobfile, "gzip " + f1)
	    logWrite(jobfile, "gzip " + f2)


    #Submit jobfile to cluster	
    os.system("chmod +x " + jobfile.name)

    if sampleNum==1:
        runAllfile = open(runAllDir + '/runAllJobs.txt',"w+");
    else:        
        runAllfile = open(runAllDir + '/runAllJobs.txt',"a+");
 
    logWrite(runAllfile, jobfile.name)    
    

    #Close files
    SampleLog.close()
    runAllfile.close()
    logfile.close()
    jobfile.close()

#Writing to logfile
def logWrite(fid, message):
    fid.write(message + '\n')


def JobFileHeader(jobfilePath):
    logWrite(jobfilePath, "#!/bin/bash -l\n");

def setDefaults(flags):
    flags['queue'] = "regevlab";
    flags['extend'] = 264;
    flags['qmem'] = 6;
    flags['qcore'] = 6;
    flags['strand'] = "NONE";
    
    procs = ['tag', 'top', 'bwt', 'pic', 'tdf', 'rsem', 'fq', 'miso', 'dup', 'trim', 'cuff'];
    for p in procs:
        flags[p + 'queue'] = flags['queue'];
	flags[p + 'mem'] = flags['qmem'];
	flags[p + 'core'] = flags['qcore'];
    
    return flags;

def runTophat(jobfile, flags, sampleName, samplesDir, bamDir, fq1=None, fq2=None):    


    #Directories
    topOutDir = samplesDir + '/' + sampleName + '/TOPHAT'
    tophatBam = topOutDir + '/' + sampleName + '.header.bam';
    finalBam =  topOutDir + '/' + sampleName + '.header.sort.bam';
    oldBam = topOutDir + '/accepted_hits.bam'
    finalBamLink = bamDir + '/' + sampleName + '.header.sort.bam';

    if not os.path.exists(topOutDir): os.makedirs(topOutDir)
    
    tophat_Path = flags['tophat_Path']
    samtools_Path = flags['samtools_Path']
       
    #Include GTF file

    if flags["single-end"] == 1:
        topArgs1 = "--num-threads %d"%(flags['topcore']);   
    else:
	topArgs1 = "--num-threads %d --mate-inner-dist %s --mate-std-dev %s"%(flags['topcore'], flags['topMateDist'], flags['topMateStd']);
    
    if 'gtf' in flags.keys():
        topArgs2 = "--G %s"%(flags['gtf'])
    if 'tophatRNA' in flags.keys():
	topArgs2 = "--transcriptome-index %s %s"%(flags['tophatTransIndex'], flags['fasta']);
    else:
	topArgs2 = flags['fasta'];

    #Tophat CMD
    tophatCMD = "%s/tophat2 -o %s %s %s %s %s"%(tophat_Path, topOutDir, topArgs1, topArgs2, fq1, fq2);
    sortCMD = "%s/samtools sort %s/accepted_hits.bam %s/%s.header.sort"%(samtools_Path, topOutDir,topOutDir, sampleName);

    bigCMD = tophatCMD + " \\\n&& " + sortCMD;
	
    #Job write
    logWrite(jobfile, '\n#Running Tophat \n');
    logWrite(jobfile, bigCMD);
    logWrite(jobfile, '\n#Linking BAM file \n');
    logWrite(jobfile,"ln -s %s %s/"%(finalBam, bamDir))
    indexCMD = "%s/samtools index %s"%(samtools_Path,finalBamLink)
    logWrite(jobfile,'\n#Indexing BAM file at link location\n');
    logWrite(jobfile,indexCMD);

    return (finalBam, finalBamLink, oldBam)


def runPicard(jobfile, flags, sampleName, metricsDir, bamName=None):
    
    picardOutDir = metricsDir + '/Picard/' + sampleName
    if not os.path.exists(picardOutDir): os.makedirs(picardOutDir)
    tempDir = picardOutDir + '/temp'
    if not os.path.exists(tempDir): os.makedirs(tempDir)
    prefix=picardOutDir + '/' + sampleName

    picard_Path = flags['picard_Path']
    
    picardCMD1_1 = "java -Xmx2g -jar %s/CollectRnaSeqMetrics.jar TMP_DIR=%s INPUT=%s OUTPUT=%s.rna_metrics.txt CHART_OUTPUT=%s.rna_coverage.pdf" % (picard_Path, tempDir, bamName, prefix, prefix)
    picardCMD1_2 = "REF_FLAT=%s STRAND=%s RIBOSOMAL_INTERVALS=%s" % (flags['refFlat'], flags['strand'], flags['riboList'])
    picardCMD1 = "%s %s" % (picardCMD1_1, picardCMD1_2)

    picardCMD2 =  "java -Xmx2g -jar %s/CollectAlignmentSummaryMetrics.jar INPUT=%s OUTPUT=%s.aln_metrics.txt" % (picard_Path, bamName, prefix)
    
    bigCMD = "%s \\\n&& %s" % (picardCMD1, picardCMD2)
    
    if not flags["single-end"]==1:
        picardCMD3 = "java -Xmx2g -jar %s/CollectInsertSizeMetrics.jar INPUT=%s OUTPUT=%s.ins_metrics.txt H=%s.ins_metrics.histogram.pdf" % (picard_Path, bamName, prefix, prefix)
	bigCMD = "%s \\\n&& %s \\\n&& %s" % (picardCMD1, picardCMD2, picardCMD3)	
      
    logWrite(jobfile, '\n#Runnning Picard tools \n')
    logWrite(jobfile, bigCMD)

def Bam2Tdf(jobfile, flags, sampleName, tdfDir, bamName=None):
    

    tdfName = tdfDir + '/' + sampleName + '.tdf'
    igvtools_Path = flags['igvtools_Path']
    
    tdfCMD = "%s/igvtools count -w 1 %s %s %s"%(igvtools_Path, bamName, tdfName, flags['igv_genome'])

    logWrite(jobfile, '\n#Converting BAM to TDF\n')
    logWrite(jobfile, tdfCMD)

def runRSEM(jobfile, flags, sampleName, samplesDir, metricsDir, RSEM_AllOutDir, fq1=None, fq2=None):

    rsemOutDir =  samplesDir + '/' + sampleName + '/RSEM'
    if not os.path.exists(rsemOutDir): os.makedirs(rsemOutDir)
    rsemSampleName = rsemOutDir + '/' + sampleName + '.rsem'
    rsemPlotOutDir = metricsDir + '/RSEM_ins_size'
    if not os.path.exists(rsemPlotOutDir): os.makedirs(rsemPlotOutDir)
    #rsem_Bam = rsemOutDir + '/' + sampleName + '.rsem.bam';
    #rsem_Sam = rsemOutDir + '/' + sampleName + '.rsem.sam';
    
    rsemPath = flags['rsem_Path']


    if flags["single-end"]==1:
        rsemArgs = "--bowtie-path %s --estimate-rspd --no-bam-output --bowtie-chunkmbs 512 -p %d %s %s %s" % (flags['bwt_Path'], flags['rsemcore'], fq1, flags['rsemTransIndex'], rsemSampleName)
    else:
        rsemArgs = "--bowtie-path %s --estimate-rspd --no-bam-output --bowtie-chunkmbs 512 -p %d --paired-end %s %s %s %s" % (flags['bwt_Path'],flags['rsemcore'], fq1,fq2, flags['rsemTransIndex'], rsemSampleName)
	    
    rsemCMD = rsemPath + "/rsem-calculate-expression " + rsemArgs
    rsemInsOut = rsemPlotOutDir  + '/' + sampleName + '_RSEM_ins_size.pdf' 
    rsemPlotCMD = rsemPath + "/rsem-plot-model %s %s" % (rsemSampleName, rsemInsOut) 
    rsemLnCMD1 = "ln -s %s/%s.rsem.genes.results %s/" % (rsemOutDir, sampleName, RSEM_AllOutDir)
    rsemLnCMD2 = "ln -s %s/%s.rsem.isoforms.results %s/" % (rsemOutDir, sampleName, RSEM_AllOutDir)

    bigCMD = "%s \\\n&& %s \\\n&& %s \\\n&& %s" % (rsemCMD, rsemPlotCMD, rsemLnCMD1, rsemLnCMD2)
    logWrite(jobfile, '\n#Runnning RSEM \n')
    logWrite(jobfile, bigCMD)

def runCufflinks(jobfile, flags, sampleName, samplesDir, Cuff_AllOutDir, bamName=None):
    
    cufflinksOutDir = samplesDir + '/' + sampleName + '/Cufflinks'
    if not os.path.exists(cufflinksOutDir): os.makedirs(cufflinksOutDir)
    
    cufflinks_Path = flags['cufflinks_Path']

    cuffArgs = "-u -p %d -g %s -b %s -L %s" % (flags['cuffcore'],flags['gtf'], flags['cuffRef'], sampleName)
    cuffCMD = "%s/cufflinks %s %s -o %s" % (cufflinks_Path, cuffArgs, bamName, cufflinksOutDir)
    cuffLnCMD1 = "ln -s %s/genes.fpkm_tracking %s/%s.cuff.genes.fpkm" % (cufflinksOutDir, Cuff_AllOutDir, sampleName)
    cuffLnCMD2 = "ln -s %s/isoforms.fpkm_tracking %s/%s.cuff.isoforms.fpkm" % (cufflinksOutDir, Cuff_AllOutDir, sampleName)
    bigCMD = "%s \\\n&& %s \\\n&& %s" % (cuffCMD, cuffLnCMD1, cuffLnCMD2)
    logWrite(jobfile, '\n#Running Cufflinks\n')
    logWrite(jobfile, bigCMD)


def runRNASeqQC(jobfile, flags, sampleName, metricsDir, bamName=None):
    
    rsqcOutDir = metricsDir +  '/RNASeqQC/' + sampleName
    if not os.path.exists(rsqcOutDir): os.makedirs(rsqcOutDir)
    
    picard_Path = flags['picard_Path']
    rsqc_Path = flags['rna_seQC_Path']

    #Add read groups
    tmpBamWrg = rsqcOutDir + '/' + os.path.basename(bamName)
    tmpBamWrg = tmpBamWrg.replace(".bam",".wrg.bam")
    rsqcCMD01 = "java -Xmx4G -jar %s/AddOrReplaceReadGroups.jar I=%s O=%s ID=none LB=none PL=none PU=none SM=none" % (picard_Path,bamName, tmpBamWrg)

    # now, reorder the file according to the reference fasta file
    tmpBamWrg_sorted = tmpBamWrg.replace(".bam", ".resorted.bam")
    rsqcCMD02 = "java -Xmx4G -jar %s/ReorderSam.jar I=%s R=%s O=%s" % (picard_Path, tmpBamWrg, flags['fasta'], tmpBamWrg_sorted)

    #Mark Duplicates
    tmpBamWrg_sorted_dups_marked = tmpBamWrg_sorted.replace(".bam",".dups_marked.bam")
    rsqcCMD03 = "java -Xmx4G -jar %s/MarkDuplicates.jar I=%s O=%s M=%s.metrics, REMOVE_DUPLICATES=false" % (picard_Path, tmpBamWrg_sorted, tmpBamWrg_sorted_dups_marked, tmpBamWrg_sorted_dups_marked)
    
    # and now run rnaseq-qc:
    # index the bam file: samtools index file.bam
    samCMD = "%s/samtools index %s" % (flags['samtools_Path'], tmpBamWrg_sorted_dups_marked)
    rsqcCMD = "java -Xmx4G -jar %s/RNA-SeQC.jar -n 1000 -s \"%s|%s|rnaseQC\" -t %s -r %s -o %s" % (rsqc_Path, sampleName, tmpBamWrg_sorted_dups_marked, flags['gtf'], flags['fasta'], rsqcOutDir)
   
    if 'rRNA_fasta' in flags.keys():
        rsqcCMD = "%s -BWArRNA %s" % (rsqcCMD, flags['rRNA_fasta'])

    bigCMD = "%s \\\n&& %s \\\n&& %s \\\n&& %s \\\n&& %s" % (rsqcCMD01, rsqcCMD02, rsqcCMD03, samCMD, rsqcCMD)
    logWrite(jobfile, '\n#Runnning RNA-seQC \n')
    logWrite(jobfile, bigCMD)

def runCountrRNA(jobfile, flags, sampleName, metricsDir, fq1=None, fq2=None):
    
    rRNA_OutDir = metricsDir + '/rRNA'
    if not os.path.exists(rRNA_OutDir): os.makedirs(rRNA_OutDir)

    bwa_Path = flags['bwa_Path']
    
    CMD1 = "%s/bwa aln %s %s > %s/%s_rRNA.left.sai" % (bwa_Path, flags['rRNA_fasta'],fq1, rRNA_OutDir, sampleName)
    if flags['single-end']==1:
        CMD2 = "%s/bwa samse %s %s/%s_rRNA.left.sai %s > %s/%s_rRNA.sam" % (bwa_Path, flags['rRNA_fasta'], rRNA_OutDir,sampleName, fq1, rRNA_OutDir, sampleName)
        CMD3 = "cat %s/%s_rRNA.sam | %s/samtools view -S -F 4 - | cut -f 1 | sort -u | wc -l > %s/%s_rRNA.sam.frag_count" % (rRNA_OutDir,sampleName, flags['samtools_Path'], rRNA_OutDir, sampleName);
	bigCMD = "%s \\\n&& %s \\\n&& %s" % (CMD1, CMD2, CMD3)
    else:
        CMD2 = "%s/bwa aln %s %s > %s/%s_rRNA.right.sai" % (bwa_Path, flags['rRNA_fasta'],fq2, rRNA_OutDir, sampleName)
        CMD3 = "%s/bwa sampe %s %s/%s_rRNA.left.sai %s/%s_rRNA.right.sai %s %s > %s/%s_rRNA.sam" % (bwa_Path, flags['rRNA_fasta'], rRNA_OutDir, sampleName, rRNA_OutDir,sampleName,fq1, fq2, rRNA_OutDir, sampleName);
	CMD4 = "cat %s/%s_rRNA.sam | %s/samtools view -S -F 4 - | cut -f 1 | sort -u | wc -l > %s/%s_rRNA.sam.frag_count" % (rRNA_OutDir,sampleName, flags['samtools_Path'], rRNA_OutDir, sampleName);
        bigCMD = "%s \\\n&& %s \\\n&& %s \\\n&& %s" % (CMD1, CMD2, CMD3, CMD4)

    logWrite(jobfile, '\n#Counting rRNA reads \n')
    logWrite(jobfile, bigCMD)

def runRSeQC(jobfile, flags, sampleName, metricsDir, bamName=None):
    
    if 'RSeQC_Path' not in flags.keys():
        print 'ERROR: Please provide RSeQC_Path'
	return

    if 'trans_bed' not in flags.keys():
        print 'ERROR: Cannot run RSeQC. transcriptome.bed missing'
        print 'Please provide a transcriptome annotation file in .bed format (Make sure it is consistent with your .gtf'
	return
	
    if 'read_len' not in flags.keys():
        print 'ERROR: Cannot run RSeQC. read length missing'
        print 'Please provide read length'
	return
    
    if 'housekeeping_bed' not in flags.keys():
        print 'WARNING: Housekeeping.bed not provided. Estimating genome coverage based on transcriptome.bed. This might take long'
	flags['housekeeping_bed']=flags['trans_bed']

    rsqc2OutDir = metricsDir + '/RSeQC'
    if not os.path.exists(rsqc2OutDir): os.makedirs(rsqc2OutDir)
    out_prefix = sampleName
    
    logWrite(jobfile, '\n#Running RSeQC \n')
    #Clipping
    Dir = rsqc2OutDir + '/clipping'
    if not os.path.exists(Dir): os.makedirs(Dir)
    cmd = "%s/clipping_profile.py -i %s -o %s/%s \n" % (flags['RSeQC_Path'], bamName, Dir, out_prefix) 
    logWrite(jobfile, cmd)

    #Mismatch
    Dir = rsqc2OutDir + '/mismatch_profile'
    if not os.path.exists(Dir): os.makedirs(Dir)
    cmd = "%s/mismatch_profile.py -i %s -l %s -o %s/%s -n 5000000 \n" % (flags['RSeQC_Path'], bamName,flags['read_len'], Dir, out_prefix) 
    logWrite(jobfile, cmd)

    #Insertion profile
    Dir = rsqc2OutDir + '/insertion_profile'
    if not os.path.exists(Dir): os.makedirs(Dir)
    cmd = "%s/insertion_profile.py -i %s -l %s -o %s/%s -n 5000000 \n" % (flags['RSeQC_Path'], bamName,flags['read_len'], Dir, out_prefix)
    logWrite(jobfile, cmd)
    
    #Deletion profile
    Dir = rsqc2OutDir + '/deletion_profile'
    if not os.path.exists(Dir): os.makedirs(Dir)
    cmd = "%s/deletion_profile.py -i %s -l %s -o %s/%s -n 5000000 \n" % (flags['RSeQC_Path'], bamName,flags['read_len'], Dir, out_prefix)
    
    logWrite(jobfile, cmd)
    
    #Inner distance
    Dir = rsqc2OutDir + '/inner_distance'
    if not os.path.exists(Dir): os.makedirs(Dir)
    cmd = "%s/inner_distance.py -i %s -o %s/%s -r %s \n" % (flags['RSeQC_Path'], bamName, Dir, out_prefix, flags['trans_bed'])
    logWrite(jobfile, cmd)

    #junction annotation
    Dir = rsqc2OutDir + '/junction_annotation'
    if not os.path.exists(Dir): os.makedirs(Dir)
    cmd = "%s/junction_annotation.py -i %s -o %s/%s -r %s \n" % (flags['RSeQC_Path'], bamName, Dir, out_prefix, flags['trans_bed'])
    logWrite(jobfile, cmd)
    
    #Junction saturation
    Dir = rsqc2OutDir + '/junction_saturation'
    if not os.path.exists(Dir): os.makedirs(Dir)
    cmd = "%s/junction_saturation.py -i %s -o %s/%s -r %s \n" % (flags['RSeQC_Path'], bamName, Dir, out_prefix, flags['trans_bed'])
    logWrite(jobfile, cmd)
    
    #Read duplication
    Dir = rsqc2OutDir + '/read_duplication'
    if not os.path.exists(Dir): os.makedirs(Dir)
    cmd = "%s/read_duplication.py -i %s -o %s/%s \n" % (flags['RSeQC_Path'], bamName, Dir, out_prefix)
    logWrite(jobfile, cmd)


    #GC content

    Dir = rsqc2OutDir + '/GC_content'
    if not os.path.exists(Dir): os.makedirs(Dir)
    cmd = "%s/read_GC.py -i %s -o %s/%s \n" % (flags['RSeQC_Path'], bamName, Dir, out_prefix)
    logWrite(jobfile, cmd)

    #Nucleotide composition bias

    Dir = rsqc2OutDir + '/NT_bias'
    if not os.path.exists(Dir): os.makedirs(Dir)
    cmd = "%s/read_NVC.py -i %s -o %s/%s\n" % (flags['RSeQC_Path'], bamName, Dir, out_prefix)
    logWrite(jobfile, cmd)

    #Read quality
    
    Dir = rsqc2OutDir + '/read_quality'
    if not os.path.exists(Dir): os.makedirs(Dir)
    cmd = "%s/read_quality.py -i %s -o %s/%s \n" % (flags['RSeQC_Path'], bamName, Dir, out_prefix)
    logWrite(jobfile, cmd)

    #RPKM saturation

    Dir = rsqc2OutDir + '/RPKM_saturation'
    if not os.path.exists(Dir): os.makedirs(Dir)
    cmd = "%s/RPKM_saturation.py -i %s -o %s/%s -r %s \n" % (flags['RSeQC_Path'], bamName, Dir, out_prefix, flags['trans_bed'])
    logWrite(jobfile, cmd)

    #Gene body coverage

    Dir = rsqc2OutDir + '/geneBody_coverage_housekeeping'
    if not os.path.exists(Dir): os.makedirs(Dir)
    cmd = "%s/geneBody_coverage.py -i %s -o %s/%s -r %s \n" % (flags['RSeQC_Path'], bamName, Dir, out_prefix, flags['housekeeping_bed'])
    logWrite(jobfile, cmd)


    #Gene body coverage

    Dir = rsqc2OutDir + '/geneBody_coverage_allgenes'
    if not os.path.exists(Dir): os.makedirs(Dir)
    cmd = "%s/geneBody_coverage.py -i %s -o %s/%s -r %s \n" % (flags['RSeQC_Path'], bamName, Dir, out_prefix, flags['trans_bed'])
    logWrite(jobfile, cmd)

def runBowtie():

    bwtArgs1 = "-q --phred33-quals -n 2 -e 99999999 -l 25 -I 1 -X 2000 -a -m 200 -S -p %i %s " % (flags['bwtcore'], flags['rsemRef']);
       
    if flags["single-end"] == 1:
        bwtCmd1 = "bowtie %s %s %s" % (bwtArgs1, f1, rsem_Sam)
	bwtCmd2 = "samtools view -bS %s > %s" % (rsem_Sam, rsem_Bam) 
    else:
	bwtCmd1 = "bowtie %s -1 %s -2 %s %s" % (bwtArgs1, f1, f2, rsem_Sam)
	bwtCmd2 = "samtools view -bS %s > %s" % (rsem_Sam, rsem_Bam)

    #write in jobfile
    logWrite(jobfile, '\n #Running Bowtie for RSEM \n')
    bigCmd = "%s \\\n&& %s \\\n&& rm %s" % (bwtCmd1, bwtCmd2, rsem_Sam) 
    logWrite(jobfile, bigCmd)

def parseInput(InputArgs, flags):

    
    #Parse Command line
    projName = InputArgs[2]
    sampleName = InputArgs[4]
    flags['single-end']= int(InputArgs[6]);
    paramsFile = InputArgs[8]
    logfile = InputArgs[10]
    

    if flags['single-end']== 1:
        f1 = InputArg[11];
	f2 = None;
    else:
        f1 = InputArgs[11];
	f2 = InputArgs[12];
    
    sampleNum = int(InputArgs[14])
    #Read Input parameters
    d = csv.reader(open(paramsFile,'r'), delimiter="\t");
    paramsList = [tuple(line) for line in d if len(line) == 2];
    for row in paramsList:
        flags[row[0]] = row[1]
    return (flags, projName, sampleName, f1, f2, logfile, sampleNum)
    
def createDirs(flags,projName):

    #Project directory
    baseBigDir = os.getcwd() + '/Results'; #This has to be consistent with main.py
    samplesDir = baseBigDir + "/samples";
    bamDir = baseBigDir + "/BAMs";
    tdfDir = baseBigDir + "/TDFs";
    cuffDir = baseBigDir + "/Cuff";
    logDir = os.getcwd() + "/logs";
    metricsDir = baseBigDir + "/metrics";
    jobDir = os.getcwd() + "/jobs";
    runAllDir = os.getcwd() + "/ClusterRun";
    RSEM_AllOutDir = os.getcwd() + "/RSEM_out";
    Cuff_AllOutDir = os.getcwd() + "/Cuff_out";

    if not os.path.exists(baseBigDir):
        os.makedirs(baseBigDir)
    if not os.path.exists(samplesDir):
        os.makedirs(samplesDir)
    if not os.path.exists(bamDir):
        os.makedirs(bamDir)
    if not os.path.exists(tdfDir):
        os.makedirs(tdfDir)
    if not os.path.exists(cuffDir):
        os.makedirs(cuffDir)
    if not os.path.exists(logDir):
        os.makedirs(logDir)
    if not os.path.exists(metricsDir):
        os.makedirs(metricsDir)	
    if not os.path.exists(jobDir):
        os.makedirs(jobDir)
    if not os.path.exists(runAllDir):
        os.makedirs(runAllDir)
    if not os.path.exists(RSEM_AllOutDir):
        os.makedirs(RSEM_AllOutDir)
    if not os.path.exists(Cuff_AllOutDir):
        os.makedirs(Cuff_AllOutDir)
    
    return (baseBigDir, samplesDir, bamDir, tdfDir, cuffDir, logDir, metricsDir, jobDir, runAllDir, RSEM_AllOutDir, Cuff_AllOutDir)
if __name__ == '__main__':
    main()
