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


# In[8]:


def SetupCheckHaplo(hvcf, reference_fasta):
    """
    If there are coordintates, use them.
    If not,they will be asked
    """

    #If you leave this empty it will be asked later
    #Coordinates:
    start = None
    end = None
    chromosome = None

    if not hvcf or reference_fasta:
        raise ValueError("You must provide the path to the vcf and the reference")
        exit(1)

    # Check if start and end are provided, if not, prompt the user
    if start is None and end is None and chromosome is None:
        start = int(input("Start: "))
        end = int(input("End: "))
        chromosome = input("Chromosome: (enter only the number)")
    else:
        raise ValueError("You must provide either all or none of the coordinates")

    return start, end, chromosome

    



# In[9]:


def DefineRange(hvcf, reference_fasta):
    
    start, end, chromosome = SetupCheckHaplo(hvcf, reference_fasta)
    print("Start: ", start)
    print("End: ", end)
    chromosome = "chr" + chromosome
    print("Chromosome: ", chromosome, "\n")

    keys = []

    with open(hvcf, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            else:
                line = line.split("\t")
                ref_start = line[1]
                #capture the element of the list starting with "END"
                ref_end = [i for i in line if i.startswith("END")]
                ref_end = ref_end[0].split("=")[1]

                if int(ref_start) <= start and int(ref_end) >= end:
                    if line[0].startswith(chromosome):
                        keys = line[4].split(",")
                        keys = [i.strip("<>") for i in keys]

    if keys == []:
        print("No range found containing the coordinates")
        return None
    else:

        print (f"there are {len(keys)} keys of ranges containing the coordinates:")
        print(keys)
        print("\n")
        return keys

                    
                    





# In[10]:


def GetCoordinates(hvcf, reference_fasta):

    keys = DefineRange(hvcf, reference_fasta)

    if keys is None:
        pass

    else:
        with open(hvcf, 'r') as f:
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


def __main__():

    """
    This script will take a vcf file and a reference fasta file and will return the coordinates of the ranges that contain the coordinates provided by the user.
    Provide the path to the vcf file and the reference fasta file.
    If executing by command line, you can provide the coordinates as arguments. If not, you will be asked for them.
    """

    import os
    import pandas as pd
    import argparse
    import sys


    if AmIaNotebook() == False:
        parser = argparse.ArgumentParser(description=__main__.__doc__)
        parser.add_argument('--hvcf', "-hf", type=str, help='The hvcf file', required=True)
        parser.add_argument('--reference_fasta', "-ref", type=str, help='The reference fasta file', required=True)
        parser.add_argument('--start', "-s", type=int, help='The start coordinate')
        parser.add_argument('--end', "-e", type=int, help='The end coordinate')
        parser.add_argument('--chromosome', "-c", type=str, help='The chromosome')
        args = parser.parse_args()

        start = args.start
        end = args.end
        chromosome = args.chromosome
        hvcf = args.hvcf
        reference_fasta = args.reference_fasta

    else:
        hvcf = "/scratch/PHG/output/vcf_files/merged_hvcfs_19092024.h.vcf"
        reference_fasta = "/scratch/PHG/data/prepared_genomes/MorexV3.fa"

    GetCoordinates(hvcf, reference_fasta)



# In[12]:


if __name__ == "__main__":
    try:
        __main__()    
    except KeyboardInterrupt:
        print("Script interrupted by user. Exiting...")
        sys.exit(0)
        raise KeyboardInterrupt

