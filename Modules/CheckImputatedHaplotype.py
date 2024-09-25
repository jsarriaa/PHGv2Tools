def CheckImputatedHaplotype():

    """
    Script to check the ammount of each pangenome genome can be found in an imputated haplotype, from an hVCF
    Prepared to work with an h.VCF from PHGv2 database. Conda can be the same as phgv2-conda.
    --help for more info
    """
    import os
    import glob
    import gzip
    import subprocess
    import matplotlib.pyplot as plt
    import argparse

    parser = argparse.ArgumentParser(description=CheckImputatedHaplotype.__doc__)
    parser.add_argument('--hvcf_folder', "-folder", type=str, help='The folder with the hvcf files')
    parser.add_argument('--hvcf_file', "-hvcf", type=str, help='The hvcf file')
    args = parser.parse_args()

    hvcf_folder = args.hvcf_folder
    hvcf_file = args.hvcf_file

    if not hvcf_folder or not hvcf_file:
        print(CheckImputatedHaplotype.__doc__)
        exit()

    def import_pangenome_genomes(hvcf_folder):

        hvcf_files = glob.glob(hvcf_folder + '/*.h.vcf.gz')

        print (f"There are {len(hvcf_files)} genomes at the pangenome")

    def list_pangenome_hvcfs(hvcf_folder):

        hvcf_files = glob.glob(hvcf_folder + '/*.h.vcf.gz')
        return hvcf_files

    def check_imputation(hvcf_files, hvcf_file):

        dict_distances = {}

        with open(hvcf_file, "r") as f:
            hvcf_lines = f.readlines()
            filtred_hvcf_lines = []
            for line in hvcf_lines:
                if "#" not in line:
                    filtred_hvcf_lines.append(line)
            
            total_ranges = len(filtred_hvcf_lines)
            print (f"Total ranges in {hvcf_file} is {total_ranges}")

        for hvcf_file in hvcf_files:
            with gzip.open(hvcf_file, "rt") as f:
                print (f"Checking {f}")
                pang_content = f.read()
                match_count = 0

                for line in filtred_hvcf_lines:
                    key = line.split("\t")[4]
                    key = key.replace("<", "")
                    key = key.replace(">", "")
                    if key == ".":
                        continue
                    if key in pang_content:
                        match_count += 1
                        #print (f"Match found for {key}")
                        continue
                    else:
                        #print (f"Match not found for {key}")
                        continue
            print (f"Match count for {hvcf_file} is {match_count} out of {total_ranges} ranges({round(match_count/total_ranges*100, 2)}%)")
            dict_distances[hvcf_file] = round(match_count/total_ranges*100, 2)

        return dict_distances

    def plot_imputation(hvcf_files, hvcf_file):

        dict_distances = check_imputation(hvcf_files, hvcf_file)

        plt.bar(range(len(dict_distances)), list(dict_distances.values()), align='center')
        plt.xticks(range(len(dict_distances)), list(dict_distances.keys()), rotation=90)
        plt.show()
        plt.savefig(hvcf_file + ".png")


    hvcf_files = list_pangenome_hvcfs(hvcf_folder)

    plot_imputation(hvcf_files, hvcf_file)



CheckImputatedHaplotype()


#TODO

# Same but based in the bp, for that:

#go for last line of haplotypes and get the last bp
#get the number of bp of each range matched, and append it to the num_of_bp_matched
#divide the num_of_bp_matched by the total bp of the last line of the haplotypes