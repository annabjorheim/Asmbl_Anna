#!/usr/bin/env python3

'''
marit.hetland@outlook.com
github marithetland
June 2019
This is an automated pipeline for use in the microbiology lab @SUS
The script takes as input short-read FASTQ files and:
1) QC's the samples with FastQC and Quast
2) Trimmes the samples with TrimGalore
3) Assembles the samples using Unicycler
4) Calculates overall coverage/sequence depth
5) Assesses species and sequence type (ST)
6) Produces a QC summary-file
'''

#import modules
import os, sys, re
import logging, time
import glob
import csv
import shutil
from shutil import copyfile
import datetime
from argparse import ArgumentParser
import pandas as pd
from pathlib import Path
from subprocess import call

#Defs
def parse_args():
    #Version
    parser = ArgumentParser(description='AMR-NGS')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + "v.1.0.0")
    
    parser.add_argument('-t','--threads', type=int, default=4, required=False, help='Specify threads to use. Default: 4')
    parser.add_argument('--noex', action='store_true', required=False, help='Do not run fastQC, multiQC, Quast, MLST or Coverage calculation.')
    parser.add_argument('--nofqc', action='store_true', required=False, help='Do not run fastQC and multiQC')
    parser.add_argument('--nomlst', action='store_true', required=False, help='Do not run MLST')
    parser.add_argument('--noquast', action='store_true', required=False, help='Do not run Quast')
    parser.add_argument('--nocov', action='store_true', required=False, help='Do not calculate X')

    parser.add_argument('--klebs', action='store_true', required=False, help='Run Kleborate')

    return parser.parse_args()


def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print ('Error: Creating directory. ' +  directory)   

def run_command(command, **kwargs):
    command_str = ''.join(command)
    #logging.info('Running shell command: {}'.format(command_str))
    try:
        exit_status = call(command_str, **kwargs)
    except OSError as e:
        message = "Command '{}' failed due to O/S error: {}".format(command_str, str(e))
        raise CommandError({"Error:": message})
    if exit_status != 0:
        message = "Command '{}' failed with non-zero exit status: {}".format(command_str, exit_status)
        raise CommandError({"Error:": message})

def file_exists(seqlist, program, path, extention):
    re_run = []
    for seq in seqlist:
        if os.path.exists(path + seq + extention):
            echo = ("File-extention " + extention + " already exists for " + seq)
        else:
            echo = ("No such file present. Running " + program + " now on " + seq)
            re_run.append(seq)
    if re_run:
        return re_run
##End defs

def main():
    
    start_time = time.time()
    args = parse_args()
    now = datetime.datetime.now()    
    todays_date = now.strftime('%Y-%m-%d_%H-%M-%S')
    
    # Set up log to stdout
    logfile= None
    logging.basicConfig(
        filename=logfile,
        level=logging.DEBUG,
        filemode='w',
        format='%(asctime)s %(message)s',
        datefmt='%m-%d-%Y %H:%M:%S')
    logging.info('Running assemblyPipe v.1.0.0')
    logging.info('command line: {0}'.format(' '.join(sys.argv)))
    
    #Arguments
    if not args.threads:
        threads=str(4)
    else:
        threads=str(args.threads)
    print('Using '+ str(threads) + ' threads')

    if args.noex:
        print('Trimming and assembling reads only, no QC or downstream analyses.')
    if args.nofqc:
        print('Trimming and assembling reads only, no QC or downstream analyses.')
    if args.nomlst:
        print('MLST will not be run.')
    if args.noquast:
        print('Quast will not be run.')
    if args.nocov:
        print('Coverage will not be calculated.')
    
    if not args.noex and not args.nofqc and not args.nomlst and not args.noquast and not args.nocov and not args.klebs:
        print("Pipeline will be run with: TrimGalore, fastQC, multiQC, Unicycler, Quast, mlst and coverage calculation.")
    if not args.noex and not args.nofqc and not args.nomlst and not args.noquast and not args.nocov and args.klebs:
        print("Pipeline will be run with: TrimGalore, fastQC, multiQC, Unicycler, Quast, mlst, coverage calculation and kleborate.")
    else:
        if args.klebs:
            print('Kleborate will be run on all samples.')

    #Set current working directory
    current_dir = os.getcwd()
    if current_dir[-1] != '/':
        current_dir = current_dir + '/'
    #Checking file extensions
    fastq_raw=(current_dir+'Fastq_raw')
    for filename in glob.glob(os.path.join(fastq_raw, '*fastq.gz')):
        shutil.move(filename, current_dir)
    if any(File.endswith(".fastq.gz") for File in os.listdir(current_dir)):
        logging.info('Input files are *.fastq.gz')
        #Rename files if named *_R?_001
        if any(File.endswith(".fastq.gz") for File in os.listdir(current_dir)):
            run_command(["rename 's/_R1_001/_1/g' ",current_dir,"*R1* && rename 's/_R2_001/_2/g' ",current_dir,"*R2*"], shell=True)
        #Add all sequences to sequence_list
        sequences = glob.glob("*fastq.gz")  #Only works for .fastq.gz suffix currently
        sequence_list = []
        for sequence in sequences:
            name = sequence.replace(".gz","")
            if name.find('_1.fastq') != -1 and name[:-8] not in sequence_list:
                sequence_list.append(name[:-8])
            if name.find('_2.fastq') != -1 and name[:-8] not in sequence_list:
                sequence_list.append(name[:-8])  

        #Check that all reads have pairs
        missing_pairs = []
        for seqName in sequence_list:
            sequence_1 = False
            sequence_2 = False
            for sequence in sequences:
                if sequence.find(seqName+'_1.f') != -1:
                    sequence_1 = True
                if sequence.find(seqName+'_2.f') != -1:
                    sequence_2 = True
            if sequence_1 == False or sequence_2 == False:
                missing_pairs.append(seqName)

        if missing_pairs != []:
            logging.info("\nNot all sequence sets have pairs:")
            for seq in missing_pairs:
                print(seq)
            logging.info("Pipeline Stopped: please fix sequence pairs\n")
            sys.exit()
        else:
            logging.info("All sequences have paired. Writing to sequence_list.txt")
            with open('sequence_list.txt', 'w') as write_sequence_list:
                for sequence in sequence_list:
                    write_sequence_list.write("%s\n" % sequence)
        
        createFolder(current_dir+'logs/') 
        #broken_files = createFolder(current_dir+'broken_files/')
        unsuccessful_sequences=[]

        createFolder(current_dir+'trimmed_reads') 
        try:
            run_command(['mv *val*gz *unpaired*gz *trimming* ./trimmed_reads  2>/dev/null'], shell=True)
        except:
            pass

        #Trimming
        run_list = []
        run_list_1 = file_exists(sequence_list, 'trimgalore', './trimmed_reads/', '_1_val_1.fq.gz')
        if run_list_1:
            for item in run_list_1:
                run_list.append(item)

        run_list_2 = (file_exists(sequence_list, 'trimgalore', './trimmed_reads/', '_2_val_2.fq.gz'))
        if run_list_2:
            for item in run_list_2:
                run_list.append(item)
                
        if run_list:
            logging.info("Running TrimGalore")
            uniq_run_list = set(run_list)
            for item in uniq_run_list: 
                try:
                    run_command(['trim_galore --paired -trim1 --retain_unpaired ',item,'_?.fastq.gz > ',current_dir,'logs/',item,'_trimgalore_',todays_date,'.log 2>&1'], shell=True)
                    logging.info(item+": TrimGalore success.")
                    run_command(['mv *trimming* *val* *unpaired* ./trimmed_reads 2>/dev/null'], shell=True)
                    #TODO:ADD size-check:run_command(['if [ -s "" ] ; then echo "WARNING: Trimmed file is empty, please check." ; fi'], shell=True)
                except:
                    logging.info(item+": Trimming unsuccessful. Removing from downstream analysis.")
                    for file in glob.glob(item+'*_fq.gz'):
                        os.remove(file)
                    sequence_list.remove(item)
                    unsuccessful_sequences.append(item)

        #Run FastQC and multiQC
        if not args.nofqc and not args.noex:
            logging.info('Running FastQC')
            createFolder(current_dir+'QC/fastQC') 
            createFolder(current_dir+'QC/multiqc_trimmed') 

            #Check if QC already exists. If it does not, run it.
            run_list = []
            fastqc_R1 = file_exists(sequence_list, 'fastqc', './QC/fastQC/', '_1_val_1_fastqc.zip')
            if fastqc_R1:
                for item in fastqc_R1:
                    run_list.append(item + ("_1_val_1.fq.gz"))
            fastqc_R2 = (file_exists(sequence_list, 'fastqc', './QC/fastQC/', '_2_val_2_fastqc.zip'))
            if fastqc_R2:
                for item in fastqc_R2:
                    run_list.append(item + ("_2_val_2.fq.gz"))
            logging.info("Running FastQC on trimmed files")
            for item in run_list: 
                try:
                    run_command(['fastqc ',current_dir,'trimmed_reads/', item, ' -o QC/fastQC > ',current_dir,'logs/',item,'_fastqc_trimmed_',todays_date,'.log 2>&1' ], shell=True)
                    logging.info(item+": FastQC success. ")
                except:
                    logging.info(item+": FastQC unsuccessful. ")
                    unsuccessful_sequences.append(item)
                
            
            #Run multiqc (will run regardless of previous versions)
            logging.info('Running MultiQC.')
            try:
                run_command(['multiqc ',current_dir,'QC/fastQC -f -o ',current_dir,'QC/multiqc_trimmed/'], shell=True)  #Add option to run only multiqc if fastqc already exists
                logging.info("MultiQC success.")
            except:
                logging.info("MultiQC failure.")
        
        #Assembly
        createFolder(current_dir+'assembly') 
        createFolder(current_dir+'assemblies/')
        createFolder(current_dir+'Fastq_raw')
        createFolder(current_dir+'success') 

        run_list = []
        for seq in sequence_list:
            if os.path.isfile(current_dir + 'success/'+seq+'_Assembly_complete.txt'):
                logging.info(seq+": Assembly complete.")
                run_command(["mv ",seq,"_?.fastq.gz Fastq_raw"], shell=True)
            else:
                logging.info(seq+": Assembly incomplete.")
                run_list.append(seq)

        trimmed_dir=(current_dir+'trimmed_reads/')
        assembly_dir=(current_dir+'assemblies/')
        
        if run_list:
            logging.info("Running Unicycler assembly on unassembled files")
            uniq_run_list = set(run_list)
            with open('uniq_run_list_as.txt', 'w') as f:
                for item in uniq_run_list:
                    f.write("%s\n" % item)
            try:
                run_command(["cd ",trimmed_dir," ; parallel --jobs ",threads," 'echo {} ; unicycler -1 {}_1_val_1.fq.gz -2 {}_2_val_2.fq.gz \
                     -o ../assembly/{}_assembly --verbosity 2 --keep 2 ; touch ../success/{}_Assembly_complete.txt; mv ../{}_?.fastq.gz ../Fastq_raw' ::: $(cat ",current_dir,"uniq_run_list_as.txt) ; cd ",current_dir], shell=True)
            except:
                logging.info(": Assembly unsuccessful.") # Removing from downstream analysis.")

        try:
            for root, dirs, files in os.walk("./assembly/"):
                if not files:
                    continue
                prefix = os.path.basename(root)
                prefix= prefix.replace('_assembly', '')
                for f in files:
                    if not f.startswith(prefix):
                        os.rename(os.path.join(root, f), os.path.join(root, "{}_{}".format(prefix, f)))
        except:
            pass
        #Copy files to assemblies-directory
        try:
            for files in glob.glob(current_dir + 'assembly/**/*assembly.fasta', recursive=True):
                filename = os.path.basename(files)
                copyfile(files, os.path.join(current_dir+'assemblies/',filename))
        except:
            pass

        #Move Reads to folders
        createFolder(current_dir+'QC/trimmed_reads')
        try:
            run_command(['mv *trimming_report.txt ./QC/trimmed_reads 2>/dev/null'], shell=True)
        except:
            pass
       
        #Run Quast
        if not args.noquast and not args.noex:
            logging.info('Running Quast on assemblies')
            createFolder(current_dir+'QC/Quast')
            try:
                run_command(['quast.py ',current_dir,'assemblies/*fasta -o ',current_dir,'QC/Quast > ',current_dir,'logs/quast_',todays_date,'.log 2>&1'], shell=True)
                logging.info("Quast successful")
                logging.info('Remember to open the transposed_report.tsv file to assess the quality of your assembled reads - main points to look at: Total contigs (<700, GC% (should match the species), total length (should match the species), and have a general look at largest contig, N50 and L50 values.')
                #Create Quast report
                quast_report = pd.read_csv(current_dir+'QC/Quast/transposed_report.tsv', sep='\t')
                for index, row in quast_report.iterrows():
                    contigs = row[1]
                    if contigs > 700:
                        print("NOTE: More than 700 contigs in "+row[0]+". Resequencing adviced.")
                    elif contigs > 400 and contigs < 700:
                        print("NOTE: More than 400 contigs in "+row[0]+". Consider resequencing.")
            except:
                logging.info("Quast unsuccessful.")
        
        #Get species and ST
        if not args.nomlst and not args.noex:
            logging.info('Looking for MLST')
            try:
                run_command(['cd ',current_dir,'assemblies/ ; mlst *fasta > mlst.tsv ; cd ',current_dir], shell= True)
                logging.info("Species and MLST identification success")
            except:
                logging.info("Species and MLST identification failed. Check input-directory. Alternatively, run mlst manually on the terminal in the ./assemblies-directory: 'mlst *fasta >> mlst.tsv '")

        #Get average coverage and std deviation
        if not args.nocov and not args.noex:
            run_list = []
            for seq in sequence_list:
                if os.path.isfile(current_dir + 'QC/Coverage/'+seq+'_Coverage.Success'):
                    logging.info(seq+": Coverage has been calculated.")
                else:
                    run_list.append(seq)
                    logging.info(seq+": Coverage has not been calculated.")
            
            if run_list:
                createFolder(current_dir+'QC/Coverage') 
                logging.info("Calculating average coverage of each sample")
                uniq_run_list = set(run_list)
                for item in uniq_run_list: 
                    logging.info(item)
                    try:
                        fasta=(current_dir+'assembly/'+item+'_assembly/'+item+'_assembly.fasta')
                        trim_1=(current_dir+'trimmed_reads/'+item+'_1_val_1.fq.gz')
                        trim_2=(current_dir +'trimmed_reads/'+item+'_2_val_2.fq.gz')
                        outfile=(current_dir+'QC/Coverage/overall_coverage.tsv') 
                        run_command(['bwa index ', fasta], shell= True)              
                        run_command(['bwa mem -t 8 ',fasta,' ',trim_1,' ',trim_2,' > input_c.sam  ; \
                            picard SamFormatConverter INPUT=input_c.sam VALIDATION_STRINGENCY=SILENT OUTPUT=input_c.bam ; \
                            picard SortSam INPUT=input_c.bam OUTPUT=input_2_c.bam VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate ; \
                            picard MarkDuplicates INPUT=input_2_c.bam VALIDATION_STRINGENCY=SILENT OUTPUT=final_cont.bam METRICS_FILE=dup_metrics ; \
                            picard BuildBamIndex INPUT=final_cont.bam VALIDATION_STRINGENCY=SILENT OUTPUT=final_cont.bam.bai ' ], shell= True)
                            
                        run_command(["echo -n '",item," \t' >> ",outfile," ; tot_size=$(samtools view -H final_cont.bam | grep -P '^@SQ' | cut -f 3 -d ':' | awk '{sum+=$1} END {print sum}') ; echo $tot_size ; samtools depth final_cont.bam | awk -v var=$tot_size '{sum+=$3; sumsq+=$3*$3} END {print sum/var \"\t\" sqrt(sumsq/var - (sum/var)**2)}' >> ",outfile," \
                            ; rm final_cont* dup_m* input* "  ], shell= True)
                        logging.info(item+": Coverage calculation success.")
                        run_command(['touch ',current_dir,'QC/Coverage/',item,'_Coverage.Success'], shell= True)
                    except:
                        logging.info(item+": Coverage-calculation unsuccessful. Removing from downstream analysis.")
                        sequence_list.remove(item)
                        unsuccessful_sequences.append(item)
            
        #Run kleborate
        #ToDO: integrate Kleborate in final report
        if args.klebs:
            try:
                run_command(['kleborate --all -a ',current_dir,'assemblies/*fasta'], shell= True) 
            except:
                print("Kleborate failed, do you have kleborate in your path?")
                pass

        logging.info("Creating lists of successful and unsuccessful sequences, see 'successful_sequences.txt' and 'failed_sequences.txt'.")
        with open("successful_sequences.txt","w") as seq_suc:
            seq_suc.write(("\n".join([str(i) for i in sequence_list] )))
        with open("failed_sequences.txt","w") as seq_unsuc:
            seq_unsuc.write(("\n".join([str(i) for i in unsuccessful_sequences] )))
        
    else:
        logging.info('ERROR: Please provide input files in fastq.gz format. If your FASTQ files are separated into one folder for each read (like when downloaded from basespace), copy and paste the following into the command line: mv */* ./ ; find . -type d -empty -delete ')

    ###Creating output-file
    seq_df = pd.DataFrame(sequence_list, columns=["Assembly"])  #Created first col with seqname for each file
    
    #need to tweak for options
    #mlst
    mlst_file = pd.read_csv(current_dir+'assemblies/mlst.tsv', sep='\t', header=None)
    mlst_file.columns = ['Assembly','species','ST','al1','al2','al3','al4','al5','al6','al7']
    mlst_df = mlst_file.replace("_assembly.fasta","", regex=True)
    mlst_df_sub = mlst_df[['Assembly','species','ST']]
    seq_df_mlst = pd.merge(seq_df, mlst_df_sub, on='Assembly', how='outer')

    quast_file = pd.read_csv(current_dir+'QC/Quast/transposed_report.tsv', sep='\t')
    quast_df = quast_file.replace("_assembly","", regex=True)
    quast_df.rename(columns={'# contigs (>= 0 bp)':'#contigs'}, inplace=True)
    quast_df.rename(columns={'Total length (>= 0 bp)':'Total_length'}, inplace=True)
    quast_df_sub = quast_df[['Assembly', '#contigs','GC (%)','N50', 'L50', 'Total_length', 'Largest contig']]
    seq_df_mlst_quast = pd.merge(seq_df_mlst, quast_df_sub, on='Assembly', how='outer')

    #AVG COV + STDEV
    cov_file = pd.read_csv(current_dir+'QC/Coverage/overall_coverage.tsv', sep='\t', header=None)
    cov_file.columns = ['Assembly','Avg_coverage','StDev']
    cov_file = cov_file.replace(" ","", regex=True)
    cov_df_sub = cov_file[['Assembly','Avg_coverage','StDev']]
    seq_df_mlst_quast_cov = pd.merge(seq_df_mlst_quast, cov_df_sub, on='Assembly', how='outer')

    #fastqc_file=
    fastqc_file = pd.read_csv(current_dir+'QC/multiqc_trimmed/multiqc_data/multiqc_fastqc.txt', sep='\t')
    fastqc_df = fastqc_file.replace("_1_val_1","", regex=True)
    fastqc_df = fastqc_df.replace("_2_val_2","", regex=True)
    fastqc_df.rename(columns={'Sample':'Assembly'}, inplace=True)
    fastqc_df.rename(columns={'Total Sequences':'#Reads'}, inplace=True)
    fastqc_df_sub = fastqc_df[['Assembly', '#Reads']]
    fastqc_df_sub=fastqc_df_sub.drop_duplicates() #All pairs should have same number of reads/sequences
    seq_df_mlst_quast_cov_fastqc = pd.merge(seq_df_mlst_quast_cov, fastqc_df_sub, on='Assembly', how='outer')

    seq_df_mlst_quast_cov_fastqc.to_csv(path_or_buf='AsmPipe_'+todays_date+'.csv', sep="\t")

    try:   
            run_command(['mv *fastq.gz ./Fastq_raw 2>/dev/null'], shell=True)
    except:
        pass
    #End of file
    total_time = time.time() - start_time
    time_mins = float(total_time) / 60
    logging.info('AMR-NGS finished in ' + str(time_mins) + ' mins.')


if __name__ == '__main__':
    main()