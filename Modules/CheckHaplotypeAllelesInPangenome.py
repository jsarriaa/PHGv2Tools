#!/usr/bin/env python
# coding: utf-8

# This script works to find to which range corresponds a sequence from the reference genome of the pangenome database built with phgv2, and check how many alleles does it have.
# 
# It needs the merged hvcf of the pangenome database, and then the coordinates

# In[7]:


def AmIaNotebook():
    try:
        shell = get_ipython().__class__.__name__
        if shell == 'ZMQInteractiveShell':
            return True  # Jupyter notebook or qtconsole
        elif shell == 'TerminalInteractiveShell':
            return False  # Terminal running IPython
        else:
            return False  # Other type of environment
    except NameError:
        return False  # Not in an interactive environment

#Check if the name of this file ends with .py or with .


def DefineRange(pangenome_hvcf, reference_fasta, start, end, chromosome):
    
    print("Start: ", start)
    print("End: ", end)
    chromosome = "chr" + chromosome
    print("Chromosome: ", chromosome, "\n")



    keys = []
    new_keys = []

    with open(pangenome_hvcf, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            else:
                line = line.split("\t")
                ref_start = line[1]
                #capture the element of the list starting with "END"
                ref_end = [i for i in line if i.startswith("END")]
                ref_end = ref_end[0].split("=")[1]

                #Check all ranges that overlap the coordinates
                if int(ref_start) >= start and int(ref_end) <= end:

                    if line[0].startswith(chromosome):
                        new_keys = line[4].split(",")
                        new_keys = [i.strip("<>") for i in new_keys]
                        keys.extend(new_keys)

                #check if a single range contains the coordinates
                elif int(ref_start) <= start and int(ref_end) >= end:
                    if line[0].startswith(chromosome):
                        new_keys = line[4].split(",")
                        new_keys = [i.strip("<>") for i in new_keys]
                        keys.extend(new_keys)


    if keys == []:
        print("No range found containing the coordinates")
        return None
    else:

        print (f"there are {len(keys)} keys of ranges containing the coordinates:")
        print(keys)
        print("\n")
        return keys

# In[10]:


def GetCoordinates(pangeonme_hvcf, reference_fasta, start, end, chromosome):

    keys = DefineRange(pangeonme_hvcf, reference_fasta, start, end, chromosome)

    if keys is None:
        pass

    else:
        with open(pangeonme_hvcf, 'r') as f:
            for line in f:
                if not line.startswith("#"):
                    continue
                else:
                    for key in keys:
                        if key in line:
                            line = line.split(",")
                            samplename = line[3]
                            samplename = samplename.split("=")[1]
                            region = line[4]
                            region = region.split("=")[1]
                            ref_range = line[6].split("=")[1]
                            ref_range = ref_range[:-1] #remove the last character ">"

                            print(f"Line: {samplename}\tRegion: {region}\tReference range:{ref_range} ")

                            

# In[11]:
if AmIaNotebook() == True:
    pangenome_hvcf = "/scratch/PHG/output/vcf_files/merged_hvcfs_19092024.h.vcf"
    reference_fasta = "/scratch/PHG/data/prepared_genomes/MorexV3.fa"
    start = None        #Add the coordinates here if you want to provide them
    end = None          #Add the coordinates here if you want to provide them
    chromosome = None   #Add the chromosome here if you want to provide it

def main(pangenome_hvcf, reference_fasta, start, end, chromosome):

    """
    This script will take a vcf file and a reference fasta file and will return the coordinates of the ranges that contain the coordinates provided by the user.
    Provide the path to the vcf file and the reference fasta file.
    If executing by command line, you can provide the coordinates as arguments. If not, you will be asked for them.
    """

    import os
    import pandas as pd
    import sys




    if start is None and end is None and chromosome is None:
        start = int(input("Start: "))
        end = int(input("End: "))
        chromosome = input("Chromosome: (enter only the number)")


    GetCoordinates(pangenome_hvcf, reference_fasta, start, end, chromosome)

    print("\nWork in progress: to ask if the user wants to download a fasta file with the sequence of the range")



# In[12]:


if __name__ == "__main__":
    try:
        main()    
    except KeyboardInterrupt:
        print("Script interrupted by user. Exiting...")
        sys.exit(0)
        raise KeyboardInterrupt

