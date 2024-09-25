def CheckHaplotypeAlleles():

    """
    This function checks the alleles of the haplotypes in the HVCFS file.
    It requires the hVCF file of the merged database, and the reference fasta file.
    Then, you will have to provide the coordinates of the region you want to check, and the chromosome.
    """

    import os
    import pandas as pd
    import argparse

    parser = argparse.ArgumentParser(description=CheckHaplotypeAlleles.__doc__)
    parser.add_argument('--hvcf', "-hf", type=str, help='The hvcf file')
    parser.add_argument('--reference_fasta', "-ref", type=str, help='The reference fasta file')
    args = parser.parse_args()

    hvcf = args.hvcf
    reference_fasta = args.reference_fasta

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


        # Check if start and end are provided, if not, prompt the user
        if start is None and end is None and chromosome is None:
            start = int(input("Start: "))
            end = int(input("End: "))
            chromosome = input("Chromosome: (enter only the number)")
        else:
            raise ValueError("You must provide either all or none of the coordinates")

        return start, end, chromosome

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

    GetCoordinates(hvcf, reference_fasta)

CheckHaplotypeAlleles()