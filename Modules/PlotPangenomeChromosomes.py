#!/usr/bin/env python
# coding: utf-8

# In[19]:


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

#Check if the name of this file ends with .py or with .ipynb


# In[20]:


def DefineHVCFs(pangenome_hvcf_folder, reference_hvcf):
    """
    Define the haplotype VCF files to be plotted
    """
    haplotype_files = [f for f in os.listdir(pangenome_hvcf_folder) if f.endswith('.h.vcf.gz')]
    if reference_hvcf in haplotype_files:
        haplotype_files.remove(reference_hvcf)

    haplotype_names = [f.split('/')[-1].split('.h.vcf.gz')[0] for f in haplotype_files]
    reference_name = reference_hvcf.split('.h.vcf.gz')[0].split('/')[-1]
    
    return haplotype_files, haplotype_names, reference_name


# In[21]:


def PangenomeColors(reference_name, haplotype_names):
    """
    Define the colors for the pangenome plot
    """

    columns_header = ("chrom", "chromStart", "chromEnd", "name", "score", "strand", "thickStart", "thickEnd", "itemRgb")

    colors = {
        "red": "255, 0, 0",
        "green": "0, 255, 0",
        "blue": "0, 0, 255",
        "yellow": "255, 255, 0",
        "cyan": "0, 255, 255",
        "magenta": "255, 0, 255",
        "brown": "165, 42, 42",
        "gray": "128, 128, 128",
        "orange": "255, 165, 0",
        "purple": "128, 0, 128",
        "pink": "255, 192, 203",
        "lime": "0, 255, 0",
        "navy": "0, 0, 128",
        "teal": "0, 128, 128",
        "olive": "128, 128, 0",
        "maroon": "128, 0, 0",
        "violet": "238, 130, 238",
        "gold": "255, 215, 0",
        "silver": "192, 192, 192",
        "beige": "245, 245, 220",
        "coral": "255, 127, 80",
        "indigo": "75, 0, 130",
        "khaki": "240, 230, 140",
        "lavender": "230, 230, 250",
    }

    ref_color = {f"{reference_name}": "0, 0, 0"}

    #Only use the first colors, equal to number of haplotype_names
    colors = dict(list(colors.items())[:len(haplotype_names)])

    # Create a new dictionary with identifiers as keys
    colors_dic = {haplotype_names[i]: color for i, color in enumerate(colors.values())}

    print(f"\nThe colors assigned to the haplotypes are: {colors_dic}\n")

    return colors_dic, ref_color, columns_header


# In[22]:


def CreateBedFile(pangenome_hvcf_folder, plots_folder, haplotype_names, reference_name, reference_hvcf, colors_dic, columns_header):
    """
    Create the bed files for each genome haplotype
    """

    for file in haplotype_names:
        #check if file already exists:
        print(f"starting to process {file}...")
        if os.path.exists(f"{pangenome_hvcf_folder}plots/{file}_haplotype.bed"):
            print(f"file {file}_haplotype.bed already exists, removing and making it again...")
            os.remove(f"{pangenome_hvcf_folder}plots/{file}_haplotype.bed")
            first_line = True
            with gzip.open(f"{pangenome_hvcf_folder}{file}.h.vcf.gz", 'rt', encoding='utf-8', errors='ignore') as vcf_file, open(f"{pangenome_hvcf_folder}plots/{file}_haplotype.bed", 'w') as output:
                lines = vcf_file.readlines()
                contig_haplotypes = [line for line in lines if not line.startswith("#")]
                for haplotype_line in contig_haplotypes:
                    haplotype_line = haplotype_line.strip().split('\t')
                    chr = haplotype_line[0]
                    start = haplotype_line[1]
                    end = haplotype_line[7].split('=')[1]
                    key = haplotype_line[4]
                    extract_key = re.search(r'<(.*?)>', key)
                    if extract_key:
                        key = extract_key.group(1)

                    if first_line == True:
                        output.write("\t".join(columns_header) + "\n")
                        first_line = False
                        output.write(f"{chr}\t{start}\t{end}\t{key}\t100\t+\t{start}\t{end}\t{colors_dic[file]}\n")
                        first_line = False
                    else:
                        output.write(f"{chr}\t{start}\t{end}\t{key}\t100\t+\t{start}\t{end}\t{colors_dic[file]}\n")
                print(f"finished processing {file}_haplotype.bed")


        else:
            os.makedirs(plots_folder, exist_ok=True)
            print(f"generating {file}_haplotype.bed...")
            first_line = True
            with gzip.open(f"{pangenome_hvcf_folder}{file}.h.vcf.gz", 'rt', encoding='utf-8', errors='ignore') as vcf_file, open(f"{pangenome_hvcf_folder}plots/{file}_haplotype.bed", 'w') as output:
                lines = vcf_file.readlines()
                contig_haplotypes = [line for line in lines if not line.startswith("#")]
                for haplotype_line in contig_haplotypes:
                    haplotype_line = haplotype_line.strip().split('\t')
                    chr = haplotype_line[0]
                    start = haplotype_line[1]
                    end = haplotype_line[7].split('=')[1]
                    key = haplotype_line[4]
                    extract_key = re.search(r'<(.*?)>', key)
                    if extract_key:
                        key = extract_key.group(1)

                    if first_line == True:
                        output.write("\t".join(columns_header) + "\n")
                        first_line = False
                        output.write(f"{chr}\t{start}\t{end}\t{key}\t100\t+\t{start}\t{end}\t{colors_dic[file]}\n")
                        first_line = False
                    else:
                        output.write(f"{chr}\t{start}\t{end}\t{key}\t100\t+\t{start}\t{end}\t{colors_dic[file]}\n")
                print(f"finished processing {file}_haplotype.bed")

    #Same for reference

    print(f"starting to process {reference_name}...")
    if os.path.exists(f"{pangenome_hvcf_folder}plots/{reference_name}_haplotype.bed"):
                print(f"file {reference_name}_haplotype.bed already exists, removing and making it again...")
                os.remove(f"{pangenome_hvcf_folder}plots/{reference_name}_haplotype.bed")
                first_line = True
    if reference_name != None:            
        with gzip.open(f"{reference_hvcf}", 'rt', encoding='utf-8', errors='ignore') as vcf_file, open(f"{plots_folder}{reference_name}_haplotype.bed", 'w') as output:
            lines = vcf_file.readlines()
            contig_haplotypes = [line for line in lines if not line.startswith("#")]
            for haplotype_line in contig_haplotypes:
                haplotype_line = haplotype_line.strip().split('\t')
                chr = haplotype_line[0]
                start = haplotype_line[1]
                end = haplotype_line[7].split('=')[1]
                key = haplotype_line[4]
                extract_key = re.search(r'<(.*?)>', key)
                if extract_key:
                    key = extract_key.group(1)

                if first_line == True:
                    output.write("\t".join(columns_header) + "\n")
                    first_line = False
                    output.write(f"{chr}\t{start}\t{end}\t{key}\t100\t+\t{start}\t{end}\t{colors_dic[file]}\n")
                    first_line = False
                else:
                    output.write(f"{chr}\t{start}\t{end}\t{key}\t100\t+\t{start}\t{end}\t{colors_dic[file]}\n")
            print(f"finished processing {file}_haplotype.bed")



                


# In[23]:


def WholeChrLenght(reference_fasta, chromosome_to_plot, reference_name, region_to_plot):

    """
    Get the whole chromosome length only if plotting the whole chromosome
    """

    #extract the list of real chromosomes
    with open(reference_fasta, 'r') as fasta_file:
        fasta_lines = fasta_file.readlines()
        chromosomes = [line for line in fasta_lines if line.startswith(">")]
        real_chromosomes = [chromosome.split(' ')[0].replace(">", "") for chromosome in chromosomes]
        print(f"real chr: {real_chromosomes}")

        chromosome_dict = {f"chr{i+1}": real_chromosomes[i] for i in range(len(real_chromosomes))}

    # Print the dictionary
    print(f"chr_dict = {chromosome_dict}")

    if region_to_plot == None:


        chromosome_found = False
        chromosome_length = 0
        print(f"The chromosome to plot is the {chromosome_to_plot}, which actually is called {chromosome_dict[chromosome_to_plot]} in {reference_name}...")

        for line in fasta_lines:
            if not line.startswith (f">{chromosome_dict[chromosome_to_plot]}") and not chromosome_found:
                continue
            if line.startswith (f">{chromosome_dict[chromosome_to_plot]}") and not chromosome_found:
                chromosome_found = True
                continue
            if chromosome_found:
                if line.startswith(">"):
                    break
                else:
                    chromosome_length = chromosome_length + len(line.strip())   
    
    else:
        chromosome_length = None

    return chromosome_length, chromosome_dict


# In[24]:


def DefineRegion(region_to_plot, chromosome_to_plot, plots_folder):
    """
    Define the region to plot
    """

    start = None
    end = None

    if region_to_plot:
        start, end = region_to_plot.split("-")
        if start > end:
            start, end = end, start
            region_to_plot = f"{start}-{end}"
        output_png = f"{plots_folder}pangenome_{chromosome_to_plot}_{region_to_plot}.png"
        region_to_plot = f"{start}-{end}"

    else:
        output_png = f"{plots_folder}pangenome_FULL_{chromosome_to_plot}.png"
    
    return start, end, output_png


# In[25]:


def PlotPangenome(plots_folder, chromosome_to_plot, chromosome_length, chromosome_dict, start, end, output_png, reference_name, haplotype_names):

    pangenome_haplotypes = [f for f in os.listdir(plots_folder) if f.endswith('_haplotype.bed')]
    if reference_name in pangenome_haplotypes:
        pangenome_haplotypes.remove(reference_name)

    print("there are this many haplotypes in the pangenome:")
    print(len(pangenome_haplotypes))
    print(pangenome_haplotypes)


    track_ini = f"{plots_folder}FULLpangenome_track.ini"

    if os.path.exists(track_ini):
        print("track.ini already exists, deleting and make a new file...")
        os.remove(track_ini)

    with open(track_ini, 'w') as track:
        for file in haplotype_names:
            track.write(f"[haplotype_{file}]\n")
            track.write(f"file = {plots_folder}{file}_haplotype.bed\n")
            print(f"loading for ploting in following order: {plots_folder}{file}_haplotype.bed")
            track.write(f"title = {file}\n")
            track.write(f"height = 0.5\n")
            track.write(f"color = bed_rgb\n")
            track.write("display = collapsed\n")
            track.write("labels = false\n")
            track.write("border_color = black\n")
            track.write("line_width = 1\n")
    #       track.write("type = bed\n")    #by default is the extension, not necessary
            track.write("\n")
            track.write("[spacer]\nheight = 0.5\n")

        track.write(f"[spacer]\nheight = 0.5\n")
        track.write(f"[haplotype_{reference_name}]\n")
        track.write(f"file = {plots_folder}{reference_name}_haplotype.bed\n")
        track.write(f"title = Reference_{reference_name}\n")
        track.write(f"height = 1\n")
        track.write(f"color = black\n")
        track.write("display = collapsed\n")
        track.write("labels = false\n")
        track.write("border_color = white\n")
        track.write("line_width = 1\n")
    #       track.write("type = bed\n")    #by default is the extension, not necessary
        track.write("\n")
        track.write("[spacer]\nheight = 0.1\n")

        track.write("[x-axis]\n")
        track.write("where = bottom\n")
        #track.write("label = true\n")
        #track.write("font_size = 10\n")
        #track.write(f"title = Mediterranean pangenome {chromosome_to_plot} using as reference MorexV3\n")
        track.write("[spacer]\nheight = 0.5\n")

    print(f"track.ini file is {track_ini}")

    region_to_plot = f"{chromosome_to_plot}:{start}-{end}" if start is not None and end is not None else None
    print(f"start: {start}, end: {end}")
    print(f"region_to_plot: {region_to_plot}")

    if region_to_plot is not None:
        print(f"plotting region {chromosome_dict[chromosome_to_plot]}:{start}-{end}")
        command_Genome_tracks = f"--tracks {track_ini} -o {output_png} --region {chromosome_dict[chromosome_to_plot]}:{start}-{end} --dpi 300 --width 75 --fontSize 15"

    else:   #by default whole chromosome
        print(f"plotting whole chromosome {chromosome_dict[chromosome_to_plot]}")
        print(f"chromosome_length: {chromosome_length}")
        command_Genome_tracks = f"--tracks {track_ini} -o {output_png} --region {chromosome_dict[chromosome_to_plot]}:1-{chromosome_length} --dpi 300 --width 75 --fontSize 15 "

    print(command_Genome_tracks)

    subprocess.run(f"pyGenomeTracks {command_Genome_tracks}", shell=True)

    print(output_png)

    print(f"output file is {output_png}")


# In[26]:


import os
import gzip
import re
import subprocess
import argparse

region_to_plot = None


if AmIaNotebook():
        
        pangenome_hvcf_folder = "/genoma/nfs/PHG/output/vcf_files/"
        chromosome_to_plot = "chr7" #give imput of which chromosome to plot
        reference_hvcf = "/genoma/nfs/PHG/vcf_dbs/hvcf_files/MorexV3.h.vcf.gz"
        reference_fasta = "/genoma/nfs/PHG/data/prepared_genomes/MorexV3.fa"
        region_to_plot = None

def main(pangenome_hvcf_folder, chromosome, reference_hvcf, reference_fasta, region):

    region_to_plot = region

    chromosome_to_plot = chromosome

    plots_folder = f"{pangenome_hvcf_folder}plots/"

    start, end, output_png = DefineRegion(region_to_plot, chromosome_to_plot, plots_folder)

    haplotype_files, haplotype_names, reference_name = DefineHVCFs(pangenome_hvcf_folder, reference_hvcf)

    colors_dic, ref_color, columns_header = PangenomeColors(reference_name, haplotype_names)

    CreateBedFile(pangenome_hvcf_folder, plots_folder, haplotype_names, reference_name, reference_hvcf, colors_dic, columns_header)

    chromosome_length, chromosome_dict = WholeChrLenght(reference_fasta, chromosome_to_plot, reference_name, region_to_plot)

    PlotPangenome(plots_folder, chromosome_to_plot, chromosome_length, chromosome_dict, start, end, output_png, reference_name, haplotype_names)


# In[27]:


if __name__ == "__main__":
    try:
        main()    
    except KeyboardInterrupt:
        print("Script interrupted by user. Exiting...")
        sys.exit(0)
        raise KeyboardInterrupt

