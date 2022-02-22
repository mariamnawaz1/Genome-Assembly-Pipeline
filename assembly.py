#!/usr/bin/env python3

import argparse
import os
import subprocess
import shlex
import shutil
import sys

parser = argparse.ArgumentParser()
parser.add_argument('data_dir',
                    help="""Directory containing the original <isolate_name>_1.fq.gz and <isolate_name>_2.fq.gz files for each isolate
                            OR, if the --skip-trimming option is given, the directory containing the trimmed and unzipped .fq files for each isolate""")
parser.add_argument('--threshold', type=int, default=5,
                    help="""Trimming quality threshold to start at. E.g. 5 starts trimming the read ends at Q5, 
                            and increases quality threshold by 1 until trimmed reads pass FASTQC""")
parser.add_argument('--skip_trimming', action='store_true', required=False,
                    help='Skip trimming. If using this option, provide the directory containing the trimmed reads as <data_dir>')
parser.add_argument('--skip_quast', action='store_true', required=False,
                    help='Skip QUAST of final assemblies.')

args = parser.parse_args()
    
def passing_fastqc(f):
    # Passes QC if Basic Statistics and Per base sequence quality are both "PASS"
    return True if [next(f).rstrip().split('\t')[0] for _ in range(2)] == ['PASS', 'PASS'] else False

def trim(data_dir):
    # Create trimmed and fastqc directories
    if os.path.exists('./trimmed'):
        sys.exit('Please change the name of <data_dir> to something other than "trimmed"')
    else:
        os.mkdir('./trimmed')
    if not os.path.exists('./fastqc'):
        os.mkdir('./fastqc')

    for root, subdirectories, files in os.walk(data_dir):
        for subdirectory in subdirectories:
            isolate_dir = os.path.join(root, subdirectory)
            isolate_name = isolate_dir.split("/")[-1]
            print(f"\n\n\nProcessing isolate {isolate_name} for trimming...\n")

            # Copy original files from data dir to trimmed dir and deflate
            trim_dir = f"./trimmed/{isolate_name}"
            shutil.copytree(isolate_dir, trim_dir)
            subprocess.call(shlex.split(f"gunzip -r {trim_dir}"))

            # Start trim quality threshold at args.threshold
            threshold = args.threshold
            while True:
                fastqc_base = f"./fastqc/{isolate_name}"
                if os.path.exists(fastqc_base):
                    shutil.rmtree(fastqc_base)
                os.mkdir(fastqc_base)

                # Call FASTQC and unzip the results
                print(f"\nRunning FASTQC on {isolate_name}...\n")
                subprocess.call(shlex.split(f"fastqc -o {fastqc_base} {trim_dir}/{isolate_name}_1.fq {trim_dir}/{isolate_name}_2.fq"))
                subprocess.call(shlex.split(f"unzip {fastqc_base}/{isolate_name}_1_fastqc.zip -d {fastqc_base}"))
                subprocess.call(shlex.split(f"unzip {fastqc_base}/{isolate_name}_2_fastqc.zip -d {fastqc_base}"))

                with open(f"{fastqc_base}/{isolate_name}_1_fastqc/summary.txt", 'r') as f1, open(f"{fastqc_base}/{isolate_name}_2_fastqc/summary.txt", 'r') as f2:
                    # Check if passing QC
                    if passing_fastqc(f1) and passing_fastqc(f2):
                        print(f"\n\n{isolate_name} is passing QC!!!\n\n")
                        break
                    else:
                        # If not passing, trim both ends with bbduk at current quality threshold
                        print(f"\n{isolate_name} not passing QC\n")
                        if os.path.exists(trim_dir):
                            shutil.rmtree(trim_dir)
                        os.mkdir(trim_dir)
                        cmd = (
                            f"bbduk.sh in1={isolate_dir}/{isolate_name}_1.fq.gz in2={isolate_dir}/{isolate_name}_2.fq.gz "
                            f"out1={trim_dir}/{isolate_name}_1.fq out2={trim_dir}/{isolate_name}_2.fq "
                            f"qtrim=rl trimq={threshold}"
                        )
                        with open(f"{trim_dir}/trim.log", 'w') as log_f:
                            print(f"\nTrimming with threshold: {threshold}...\n")
                            # Call bbduk and output trimming logs to trim.log
                            subprocess.call(shlex.split(cmd), stdout=log_f, stderr=log_f)
                            print(f"\nDone trimming with threshold: {threshold}\n")
                        
                        threshold += 1

    print("Generating multiqc report...")
    subprocess.call(shlex.split('multiqc ./fastqc'))
    print("\n\n\nDone trimming!\n\n\n")

def create_assemblies(trimmed_dir):
    if not os.path.exists('./spades_outputs'):
        os.mkdir('./spades_outputs')

    for root, subdirectories, files in os.walk(trimmed_dir):
        for subdirectory in subdirectories:
            isolate_dir = os.path.join(root, subdirectory)
            isolate_name = isolate_dir.split("/")[-1]
            print(f"\nAssembling isolate {isolate_name} with SPAdes\n")
            
            trim_base = f"{trimmed_dir}/{isolate_name}"

            # Run SPAdes assembly
            subprocess.call(shlex.split(f"spades.py --pe1-1 {trim_base}/{isolate_name}_1.fq --pe1-2 {trim_base}/{isolate_name}_2.fq -o ./spades_outputs/{isolate_name}/"))
            print(f"\nSPAdes assembly for {isolate_name} complete\n")
            
def copy_final_assemblies(spades_output_dir):
    if not os.path.exists('./assemblies'):
        os.mkdir('./assemblies')
 
    for root, subdirectories, files in os.walk(spades_output_dir):
        for subdirectory in subdirectories:
            if 'CGT' in subdirectory:
                isolate_dir = os.path.join(root, subdirectory)
                isolate_name = isolate_dir.split("/")[-1]
 
                print(f"\nCopying contigs.fasta for {isolate_name} to assemblies directory\n")
                copy_path = f"./assemblies/{isolate_name}"
                if not os.path.exists(copy_path):
                    os.mkdir(copy_path)
                shutil.copy(f"{spades_output_dir}/{isolate_name}/contigs.fasta", copy_path)

def quast(assemblies_root):
    for root, subdirectories, files in os.walk(assemblies_root):
        for subdirectory in subdirectories:
            if 'CGT' in subdirectory:
                isolate_dir = os.path.join(root, subdirectory)
                isolate_name = isolate_dir.split("/")[-1]

                print(f"\nRunning QUAST for {isolate_name}\n")
                isolate_base = f"{assemblies_root}/{isolate_name}"
                subprocess.call(shlex.split(f"quast -o {isolate_base}/quast_out {isolate_base}/contigs.fasta"))


if __name__ == "__main__":
    if not args.skip_trimming: # Trim reads before assembling
        trim(args.data_dir)
        trimmed_dir = './trimmed'
    else: # Skip trimming
        trimmed_dir = args.data_dir.split('/')[0]

    # Create assemblies
    create_assemblies(trimmed_dir)

    # Copy final assemblies to ./assemblies directory
    copy_final_assemblies('./spades_outputs')

    # Run QUAST on final assemblies
    if not args.skip_quast:
        quast('./assemblies')

    print("\n\n\nDONE - EXITING\n\n\n")
