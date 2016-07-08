'''
Created February 5, 2016

@author: sbrown

This script ensures that every read to be used as input for MiTCR is only comprised of ACTGN, and ensures reads are at least 40 bases long.
Input: Path to directory containing fastq files to process.

Python requirements: v2.x.x

'''
import re
import sys
#import subprocess
#import time
import os

## requires python 2.x.x


i = 0
block = ""
valid = True
reg = re.compile('^[ACTGN]+$')


numReads = 0
d = sys.argv[1]
outDir = sys.argv[2]


reads = open(os.path.join(outDir,"reads.fq"),"w")



for (root,dir,files) in os.walk(d):
    for f in files:
        if (f.endswith(".fastq") or f.endswith(".fq")) and f != "reads.fq":
            for line in open(f, "r"):
                i += 1
                block += line.rstrip() + "\t"
                if i == 2:
                    if(reg.match(line.upper()) and len(line)>40):
                        valid = True
                    else:
                        valid = False
                elif i == 4 and valid:
                    reads.write(block.replace("\t","\n"))   # cleaves off trailing \t
                    block = ""
                    i = 0
                    numReads += 1
                elif i == 4 and not valid:
                    #outfile.write(block)
                    block = ""
                    i = 0
            
            
            ## For initial application, this was run on copied fastq files.
            ## Space was an issue, so after a fastq had been processed, and
            ## valid reads written to a new file, the fastq was deleted.
            ## Removed this functionality from this version of script.
            '''
            VALID_COM = False
            while not VALID_COM:
                cmd = "rm " + f
                call = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                (res, err) = call.communicate()
                if err.decode('ascii') == "":
                    VALID_COM = True
                else:
                    print("ERROR DELETING " + f + ".\n")
                    #logMessage("Waiting 1 minute to try again...\n",logFileName)
                    time.sleep(60)
             '''       
    
    
outfile = open(os.path.join(outDir,"numReads.txt"),"w")
outfile.write(str(numReads)+"\n")
outfile.close()
    
    
reads.close()
