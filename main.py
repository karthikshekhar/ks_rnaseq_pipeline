import csv
import sys
import os

#os.system('./bash_header.sh')
gridRunnerPath = "/seq/regev_genome_portal/SOFTWARE/HpcGridRunner/hpc_cmds_GridRunner.pl"


#Project Name
projName = 'habenula1'
headDir = 'habenula'
logfile = projName + '_logfile.txt'
clusterRunFilePath = os.getcwd() + '/' + headDir + '/' +  projName + "/ClusterRun/runAllJobs.txt";
gridConfPath =   os.getcwd() + '/' + headDir + '/' +  projName + "/hpc_run.conf";

#Get samples file
samplesFile = 'samples.txt'

#parse samples file into samples
d = csv.reader(open(samplesFile,'r'),delimiter="\t");
fileList = [tuple(line) for line in d];

if len(fileList[1])==2:
    singleEndFlag = 1;
    pairedEndFlag = 0;
else:
    singleEndFlag = 0
    pairedEndFlag = 1;

l=0
for f in fileList:
    
    
    if pairedEndFlag == 1:
        f1=f[1]; f2=f[2]
    	fName = f[0]
        #Call pipeline
        l=l+1;
        cmd = "python example.py --project " + projName + " --sample " + fName + " --single-end " + str(singleEndFlag) + " --paramsFile params.txt --logFile " + logfile + " " + f1 + " " + f2 + " --num " + str(l)
        os.system(cmd)
    
    if singleEndFlag == 1:
        f1=f[1]; fName = f[0]
	#Call pipeline
	l=l+1;
	cmd = "python example.py --project " + projName + " --sample " + fName + " --single-end " + str(singleEndFlag) + " --paramsFile params.txt --logFile " + logfile + " " + f1 + " " + None + " --num " + str(l)
        os.system(cmd)   
    

#Collect QC stats + expression values from all files
summDir = os.getcwd() + '/Summary'
if not os.path.exists(summDir): os.makedirs(summDir)

os.system('echo -e "#!/bin/bash \n" > genSummary.sh')
os.system('echo -e "\nQC stats\n" >> genSummary.sh')
shortQC_CMD = "/ahg/regevdata/users/karthik/SOFTWARE/ks_rnaseq_pipeline/QC_summary/short_summary_stats.pl --reads_list_file %s --project_base_dir %s > %s/QCshort.txt" % (samplesFile,os.getcwd(), summDir);
os.system('echo -e "%s\n" >> genSummary.sh'%(shortQC_CMD))
longQC_CMD = "/ahg/regevdata/users/karthik/SOFTWARE/ks_rnaseq_pipeline/QC_summary/long_summary_stats.pl %s %s > %s/QClong.txt" % (samplesFile,os.getcwd(), summDir);
os.system('echo -e "%s\n" >> genSummary.sh'%(longQC_CMD))

os.system('echo -e "\nGet genes and isoforms list for RSEM" >> genSummary.sh')
CMDgene = "find %s/RSEM_out/ | grep genes.results > %s/rsem.genes.list" % (os.getcwd(),summDir)
CMDisoforms="find %s/RSEM_out/ | grep isoforms.results > %s/rsem.isoforms.list" % (os.getcwd(),summDir)
os.system('echo -e "%s" >> genSummary.sh'%(CMDgene))
os.system('echo -e "%s" >> genSummary.sh'%(CMDisoforms))

os.system('echo -e "\nFPKM, TPM and counts matrix for RSEM\n" >> genSummary.sh')
CMDfpkm_gene = "/ahg/regevdata/users/karthik/SOFTWARE/ks_rnaseq_pipeline/exp_quantification/merge_RSEM_output_to_matrix.pl --rsem_files %s/rsem.genes.list --mode fpkm > %s/rsem.genes.fpkm.matrix" % (summDir, summDir) 
CMDfpkm_isoform = "/ahg/regevdata/users/karthik/SOFTWARE/ks_rnaseq_pipeline/exp_quantification/merge_RSEM_output_to_matrix.pl --rsem_files %s/rsem.isoforms.list --mode fpkm > %s/rsem.isoforms.fpkm.matrix"  % (summDir, summDir) 
CMDtpm_gene = "/ahg/regevdata/users/karthik/SOFTWARE/ks_rnaseq_pipeline/exp_quantification/merge_RSEM_output_to_matrix.pl --rsem_files %s/rsem.genes.list --mode tpm > %s/rsem.genes.tpm.matrix"  % (summDir, summDir) 
CMDtpm_isoform = "/ahg/regevdata/users/karthik/SOFTWARE/ks_rnaseq_pipeline/exp_quantification/merge_RSEM_output_to_matrix.pl --rsem_files %s/rsem.isoforms.list --mode tpm > %s/rsem.isoforms.tpm.matrix" % (summDir, summDir) 
CMDcounts_gene = "/ahg/regevdata/users/karthik/SOFTWARE/ks_rnaseq_pipeline/exp_quantification/merge_RSEM_output_to_matrix.pl --rsem_files %s/rsem.genes.list --mode counts > %s/rsem.genes.counts.matrix" % (summDir, summDir) 
CMDcounts_isoform = "/ahg/regevdata/users/karthik/SOFTWARE/ks_rnaseq_pipeline/exp_quantification/merge_RSEM_output_to_matrix.pl --rsem_files %s/rsem.isoforms.list --mode counts > %s/rsem.isoforms.counts.matrix" % (summDir, summDir) 
bigCMD = "%s \n %s \n %s \n %s \n %s \n %s" % (CMDfpkm_gene, CMDfpkm_isoform, CMDtpm_gene, CMDtpm_isoform, CMDcounts_gene, CMDcounts_isoform)
os.system('echo -e "%s\n" >> genSummary.sh'%(bigCMD))

os.system('echo -e "\nGet genes and isoforms list for Cufflinks\n" >> genSummary.sh')
CMDgene = "find %s/Cuff_out/ | grep genes.fpkm > %s/cuff.genes.list" % (os.getcwd(), summDir)
os.system('echo -e "%s" >> genSummary.sh'%(CMDgene))
CMDisoforms = "find %s/Cuff_out/ | grep isoforms.fpkm > %s/cuff.isoforms.list" % (os.getcwd(), summDir)
os.system('echo -e "%s" >> genSummary.sh'%(CMDisoforms))
os.system('echo -e "\nFPKM matrix for Cufflinks\n" >> genSummary.sh')
CMDfpkm_gene = "/ahg/regevdata/users/karthik/SOFTWARE/ks_rnaseq_pipeline/exp_quantification/merge_cuff_fpkm_to_matrix.pl %s/cuff.genes.list > %s/cuff.genes.fpkm.matrix" % (summDir, summDir) 
CMDfpkm_isoform = "/ahg/regevdata/users/karthik/SOFTWARE/ks_rnaseq_pipeline/exp_quantification/merge_cuff_fpkm_to_matrix.pl %s/cuff.isoforms.list > %s/cuff.isoforms.fpkm.matrix" % (summDir, summDir) 
bigCMD = "%s \n %s" % (CMDfpkm_gene, CMDfpkm_isoform)
os.system('echo -e "%s\n" >> genSummary.sh'%(bigCMD))

os.system("chmod +x genSummary.sh")

#Run on grid
#gridRunCMD = gridRunnerPath + " -c " + clusterRunFilePath + " -G " + gridConfPath 
#os.system(gridRunCMD)
