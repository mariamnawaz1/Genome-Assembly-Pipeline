#!/usr/bin/env python3

import os
import subprocess
import shlex
import shutil

def passing_fastqc(f):
    # Passes QC if Basic Statistics and Per base sequence quality are both "PASS"
    return True if [next(f).rstrip().split('\t')[0] for _ in range(2)] == ['PASS', 'PASS'] else False

# Create trimmed and fastqc directories
if not os.path.exists('./trimmed'):
    os.mkdir('./trimmed')
if not os.path.exists('./fastqc'):
    os.mkdir('./fastqc')

data_dir = '/home/groupa/trimming/data'
for root, subdirectories, files in os.walk(data_dir):
    for subdirectory in subdirectories:
        isolate_dir = os.path.join(root, subdirectory)
        isolate_name = isolate_dir.split("/")[-1]
        print(f"\n\n\nProcessing isolate {isolate_name}...\n")

        # Copy original files from data dir to trimmed dir and deflate
        trim_dir = f"./trimmed/{isolate_name}"
        shutil.copytree(isolate_dir, trim_dir)
        subprocess.call(shlex.split(f"gunzip -r {trim_dir}"))

        # Start trim quality threshold at Q5
        threshold = 5
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


print("\n\n\nDone trimming!\n\n\n")

print("Generating multiqc report...")
subprocess.call(shlex.split('multiqc ./fastqc'))
print("Done")

