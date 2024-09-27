#!/usr/bin/env python
# coding: utf-8

# Code to process an h.vcf imputated at the database of a Practical haplotype Graph v2 to plot the ideogram
# 
# Performed from output of:
# phg version 2.4.4.158
# 
# This code has been developed at 19 sept 2024
# Is not in the github repository
# 
# run with pygenometracks conda, in python env
# 

# In[61]:


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


# In[62]:


def assign_colors(genomes, reference_name):
    """
    Assign a color to each genome in the pangenome, up to 20 genomes at this moment. Reference is also assigned as Black
    You need to provide the list of genomes and the reference genome name    
    """
    colors = {
    "red": "255,0,0",
    "green": "0,255,0",
    "blue": "0,0,255",
    "yellow": "255,255,0",
    "cyan": "0,255,255",
    "magenta": "255,0,255",
    "brown": "165,42,42",
    "gray": "128,128,128",
    "orange": "255,165,0",
    "purple": "128,0,128",
    "pink": "255,192,203",
    "lime": "0,255,0",
    "navy": "0,0,128",
    "teal": "0,128,128",
    "olive": "128,128,0",
    "maroon": "128,0,0",
    "violet": "238,130,238",
    "gold": "255,215,0",
    "silver": "192,192,192",
    "beige": "245,245,220",
    "coral": "255,127,80",
    "indigo": "75,0,130",
    "khaki": "240,230,140",
    "lavender": "230,230,250",
    }


    colors = dict(list(colors.items())[:len(genomes)])
    colors_dic = {genomes[i]: color for i, color in enumerate(colors.values())}

    if reference_name not in genomes:
        ref_color = {f"{reference_name}": "0,0,0"}

    print(f"Reference genome color: {ref_color}, saved as ref_color")
    print(f"Genomes colors: {colors_dic}, saved as colors_dic")

    return colors_dic, ref_color


#colors_dict, ref_color = assign_colors(genomes, reference_name)

#print (f"Colors assigned to genomes: {colors_dict}")



# In[63]:


def hvcf2bed(hvcf, plot_folder, colors_dic):
    """
    This function will convert the hvcf file to a bed file, and store it at the plot folder. BED file will have the same name as hvcf. 
    It will use the coordinates of the reference.
    
    The columns are:
    1. NameOfHVCF       #Doing this for ploting software, is not strictly accurate with BED format
    2. Start
    3. End
    4. Chr             #Doing this for ploting and filtrering in this software
    5. Score (actually a dummy score)
    6. Strand
    7. thickStart
    8. thickEnd
    9. itemRgb

    *Notice that the new files WILL NOT CONTAIN in each line to which genome is the range belonging to.
    So in further steps, you must convert from the color dictionary the RGB color to each genome of the pangenome :)
    """

    if not os.path.exists(plot_folder):
        raise Exception(f"Folder {plot_folder} does not exist")
    
    hvcf_BED = (f"{plot_folder}TEMP_{os.path.basename(hvcf).split('.')[0]}.bed")

    if os.path.exists(hvcf_BED):
        #print(f"File {hvcf_BED} already exists, removing it and building it again") 
        os.remove(hvcf_BED)

    print(f"Converting {hvcf} to BED file")

    if hvcf.endswith('.gz'):                        ##Checks if the file is compressed or not and acts depending on it
        open_func = lambda f: gzip.open(f, 'rt')
    else:
        open_func = open

    with open_func(hvcf) as hvcf_file, open(hvcf_BED, 'w') as bed_file:
        first_line = True
        lines = hvcf_file.readlines()
        for line in lines:
            
            if line.startswith ('##ALT=<ID'):
                line = line.strip().split(',')
                range = line[-1]                   #take last element. Sometimes there is more than one range at the haplotype itself, so we take the reference one
                                                #Actually this bug should be fixed at the PHG version 2.4.8.162, but lets keep this code for older hvcf files
                range = range.split('=')[1]
                chr = range.split(':')[0]
                chr = chr.split ("_")[0]
                start = range.split(':')[1].split('-')[0]
                end = range.split(':')[1].split('-')[1]
                end = end.split('>')[0]
                name = line[1].split(':')[1]
                name = name.split('"')[0].replace(" ", "")
                #print(name, chr, start, end, name)

                if first_line == True:
                    bed_file.write("chrom\tchromStart\tchromEnd\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\n")
                    first_line = False
                    bed_file.write(f"{os.path.basename(hvcf).split('.')[0]}\t{start}\t{end}\t{chr}\t100\t+\t{start}\t{end}\t{colors_dic[name]}\n")
                    #print(line)
                else:
                    bed_file.write(f"{os.path.basename(hvcf).split('.')[0]}\t{start}\t{end}\t{chr}\t100\t+\t{start}\t{end}\t{colors_dic[name]}\n")
                    #print(chr)
                    pass
            

        print(f"finished processing {hvcf} to BED file {hvcf_BED}")
      



# In[64]:


def CaptureChrFromBED (hvcf, plot_folder):
    """
    This function will capture the chromosome from the BED file, storing them in a list
    """
    hvcf_BED = (f"{plot_folder}TEMP_{os.path.basename(hvcf).split('.')[0]}.bed")
    chr_list = []
    with open (hvcf_BED, 'r') as f:
        lines = f.readlines()
        lines = lines[1:]
        for line in lines:
            chr = line.split('\t')[3]
            if chr not in chr_list:
                chr_list.append(chr)
    chr_list.sort()
    print(f"Chromosomes found in the BED file are: {chr_list}")
    return chr_list


# In[65]:


def BEDSplitByChr (hvcf, plot_folder):
    pass
    """
    splits and sorts the BED file in different files by chromosome, sorting it by start position. 
    This is only in order to plot it with pygenometracks. Takes by output the BED file generated by hvcf2bed.

    """
    hvcf_BED = (f"{plot_folder}TEMP_{os.path.basename(hvcf).split('.')[0]}.bed")

    chr_list = CaptureChrFromBED(hvcf, plot_folder)
    print("Command CaptureChrFromBED executed")

    for chr in chr_list:    #check if the files splited by chr already exist, and if they do, remove them            
        chr_file = (f"{plot_folder}TEMP_{os.path.basename(hvcf).split('.')[0]}_{chr}.bed")
        if os.path.exists(chr_file):
            print(f"File {chr_file} already exists, removing it and building it again") 
            os.remove (chr_file)

    try:
        open (hvcf_BED, 'r')
    except FileNotFoundError:
        print(f"File {hvcf_BED} not found")
 
    with open (hvcf_BED, 'r') as f:
        lines = f.readlines()
        if len(lines) == 0:
            raise Exception(f"File {hvcf_BED} is empty")
        if len(chr_list) == 0:
            raise Exception(f"Chromosome list is empty")
        
        print(f"opening {hvcf_BED} to split it by chromosome")
        
        for chr in chr_list:
            chr_file = (f"{plot_folder}TEMP_{os.path.basename(hvcf).split('.')[0]}_{chr}.bed")
            with open (chr_file, 'w') as chr_bed:
                chr_bed.write("chrom\tchromStart\tchromEnd\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\n")
                for line in lines:
                    if chr in line:
                        chr_bed.write(line)
            try:
                open(chr_file, 'r')
            except FileNotFoundError:
                print(f"File {chr_file} not found") 
            else:
                if len(open(chr_file).readlines()) <= 1:  
                    raise Exception(f"File {chr_file} of the {chr} is empty\n Im actually trying to find the chromosome {chr} in the file {hvcf_BED}")
                print(f"Finished splitting {hvcf_BED} by chromosome to {chr}")


# In[66]:


def sortBED (hvcf, plot_folder):

    """
Function to sort the files, inside hvcf plotting.
Requires the hvcf file, ploting folder, and calls the function CaptureChrFromBed
    """

    chr_list = CaptureChrFromBED(hvcf, plot_folder)


#Now sort the files by start position
    for chr in chr_list:
        chr_file = (f"{plot_folder}TEMP_{os.path.basename(hvcf).split('.')[0]}_{chr}.bed")
        sorted_chr_file = (f"{plot_folder}TEMP_{os.path.basename(hvcf).split('.')[0]}_{chr}_sorted.bed")
        if os.path.exists(sorted_chr_file):
            print(f"File {sorted_chr_file} already exists, removing it and building it again") 
            os.remove(sorted_chr_file)

        try:    
            open (chr_file, 'r')
        except FileNotFoundError:
            print(f"File {chr_file} not found") 
        else:
            lines = len(open(chr_file).readlines())
            if lines == 0:
                raise Exception(f"File {chr_file} is empty")

        with open (chr_file, 'r') as f:
            lines = f.readlines()
            lines = lines[1:]
            lines.sort(key=lambda x: int(x.split('\t')[1]))
            with open (sorted_chr_file, 'w') as sorted_chr_bed:
                sorted_chr_bed.write("chrom\tchromStart\tchromEnd\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\n")
                for line in lines:
                    sorted_chr_bed.write(line)

                    try:
                        line = line.split('\t')
                        len(line) == 9  
                    except:
                        raise Exception(f"Line {line} does not have 9 columns")

            try:
                with open(sorted_chr_file, 'r') as f:
                    lines = f.readlines
            except FileNotFoundError:
                print(f"File {sorted_chr_file} not found")
            else:
                if len(open(sorted_chr_file).readlines()) <= 1:
                    raise Exception(f"File {sorted_chr_file} is empty")


        try:
            with open(sorted_chr_file, 'r')  as f:
                lines = f.readlines
        except FileNotFoundError:
            print(f"File {sorted_chr_file} not found")   
        else:
            lines = len(open(sorted_chr_file).readlines())
            if lines == 0:
                raise Exception(f"File {sorted_chr_file} is empty")
        #print(f"Finished sorting {chr_file} to {sorted_chr_file}, and removing the first one")


# In[67]:


def PrepareTrack_Hvcf2plot(hvcf, reference_name, plot_folder):
    """
    This function will prepare the tracks for pygenometracks, so it can plot the hvcf files. Indicate the hvcf initial file, and the reference name, as long
    as the folder to save the plot.
    However, it actually uses the bed files from BEDSplitByChr function.

    """

    chr_list = CaptureChrFromBED(hvcf, plot_folder)
    print(chr_list)

    #Now we prepare the tracks for pygenometracks

    track_ini = (f"{plot_folder}TEMP_{os.path.basename(hvcf).split('.')[0]}_track.ini")

    if os.path.exists(track_ini):
        print(f"File {track_ini} already exists, removing it and building it again") 
        os.remove(track_ini)

    with open (track_ini, 'w') as track:
        
        for chr in chr_list:
            track.write(f"[haplotype_{chr}]\n")
            track.write(f"title = {chr}\n")
            track.write(f"file = {plot_folder}TEMP_{os.path.basename(hvcf).split('.')[0]}_{chr}_sorted.bed \n")
            track.write(f"height = 3\n")
            track.write(f"color = bed_rgb\n")
            track.write("display = collapsed\n")
            track.write("labels = true\n")
            track.write("label_fontsize = 10\n")
            track.write("border_color = black\n")
            track.write("line_width = 0.05\n")
            track.write("\n")
            track.write("[spacer]\nheight = 0.5\n")

        track.write("[x-axis]\n")
        track.write("where = bottom\n")
        #track.write("label = true\n")
        #track.write("font_size = 10\n")
        #track.write(f"title = Mediterranean pangenome {chromosome_to_plot} using as reference MorexV3\n")
        track.write("[spacer]\nheight = 0.5\n")

    print(f"Finished preparing the tracks for pygenometracks at {track_ini}")
    return track_ini




# In[68]:


def SetExtremesToPlot(hvcf, plot_folder):
    """
    """    
    #capturing the extreme values to plot of each chromosome
    chr_list = CaptureChrFromBED(hvcf, plot_folder)
    print(chr_list)
    chr_max = None
    for chr in chr_list:
        chr_file = (f"{plot_folder}TEMP_{os.path.basename(hvcf).split('.')[0]}_{chr}_sorted.bed")
        with open (chr_file, 'r') as f:
            #print (f"Processing {chr_file}")
            lines = f.readlines()
            if lines == []:
                raise Exception(f"Empty file {chr_file}")
                
            try:    
                last_line = lines[-1]
            except:
                print(f"This is the file im working now: {lines}\n")
                raise IndexError(f"I can not find the last line in {chr_file}")

            last_line = lines[-1]
            last_coord = last_line.split('\t')[2]
            if chr_max == None or chr_max < last_coord:
                chr_max = last_coord
    chr_max = int(int(chr_max)*1.005)
    print(f"Max coord to plot will be {chr_max}")

    return chr_max


# In[69]:


def RunPygenometracksCommandHVCF(hvcf, reference_name, plot_folder, genomes):

    """
    Building and executing the command for plotting:
    It calls the functions SetExtremesToPlot and PrepareTrack_HvcfPlot
    """

    chr_max = SetExtremesToPlot(hvcf, plot_folder)
    track_ini = PrepareTrack_Hvcf2plot(hvcf, reference_name, plot_folder)
    chr_list = CaptureChrFromBED(hvcf, plot_folder)

    width = len(genomes)*12
    height = len(chr_list)*10

    command = (f"pyGenomeTracks --tracks {track_ini} --region {os.path.basename(hvcf).split('.')[0]}:1-{chr_max} --outFileName {plot_folder}/{os.path.basename(hvcf)}.png --width {width} --dpi 500")

    try: 
        result = subprocess.run(command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        print(f"running the command {command}")
        #print(result.stdout)
        #print(result.stderr)
        if "No valid intervals were found" in result.stdout or "No valid intervals were found" in result.stderr:
            raise Exception(f"Command {command} failed with error: {result.stderr}")
    except subprocess.CalledProcessError as e:
        raise Exception(f"Command {command} failed with error: {e.stderr}")
        
    print(f"the plot is saved at {plot_folder}{os.path.basename(hvcf).split('.')[0]}.png")


# In[70]:


def AddLegendToPygenometracks (hvcf, plot_folder, genomes, reference_name):
    """
    This function will add a legend to the pygenometracks plot, with the colors assigned to each genome.
    """

    print ("lets add a legend to the plot")

    # Transform the colors RGB
    colors_dict, ref_color = assign_colors(genomes, reference_name)

    # Convert colors from RGB to hex
    for genome in colors_dict:
        r, g, b = map(int, colors_dict[genome].split(','))
        colors_dict[genome] = f'#{r:02X}{g:02X}{b:02X}'
    r, g, b = map(int, ref_color[reference_name].split(','))
    ref_color[reference_name] = f'#{r:02X}{g:02X}{b:02X}'

    # Create legend elements
    legend_elements = [Line2D([0], [0], marker='o', color='w', label=genome, 
                            markerfacecolor=colors_dict[genome], markersize=100) 
                    for genome in colors_dict]
    legend_elements.append(Line2D([0], [0], marker='o', color='w', label=reference_name, 
                                markerfacecolor=ref_color[reference_name], markersize=100))

    # Create a new figure for the legend
    fig, ax = plt.subplots(figsize=(10, 10))  # Adjust size as necessary
    ax.legend(handles=legend_elements, loc='center', frameon=False,fontsize=140)

    # Remove axes for a clean legend
    ax.axis('off')

    # Save the legend as an image
    plt.savefig(f'{plot_folder}TEMP_legend.png', bbox_inches='tight', transparent=True)
    plt.close()

    # Open the main plot and the legend, but only if you are in a notebook
    plot_image = Image.open(f'{plot_folder}{os.path.basename(hvcf)}.png')
    legend_image = Image.open(f'{plot_folder}TEMP_legend.png')
    if AmIaNotebook():
        display(plot_image)
        display(legend_image)

    # Get the dimensions of the plot and legend
    plot_width, plot_height = plot_image.size
    legend_width, legend_height = legend_image.size

    # Define the position where you want to place the legend (e.g., bottom-right corner)
    position = (plot_width - legend_width - 10, plot_height - legend_height - 10)  # 10-pixel padding

    # Paste the legend onto the plot
    plot_image.paste(legend_image, position, legend_image)

    # Save the final merged image
    plot_image.save(f'{plot_folder}{os.path.basename(hvcf)}.png')

    print ("Done! The legend is added to the plot\n")
    print(f"Your plot is saved at {plot_folder}{os.path.basename(hvcf)}.png")


# In[71]:


import os
import gzip
import re
import subprocess
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from PIL import Image
from matplotlib.lines import Line2D
import argparse
import sys

def main():
    """
    This function will plot the hvcf file using pygenometracks. It will use the track_ini file generated by PrepareTrack_Hvcf2plot
    """
    if AmIaNotebook() == False:
        parser = argparse.ArgumentParser(description=main.__doc__)
        parser.add_argument("--input-hvcf", "-hvcf", help="Path to the hvcf file")
        parser.add_argument("--pangenome-hvcf-folder", "-folder", help="Path to the folder where the plots will be saved")
        parser.add_argument("--reference-hvcf", "-ref", help="Path to the reference hvcf file")
        args = parser.parse_args()

        hvcf = args.input_hvcf
        pangenome_hvcf = args.pangenome_hvcf_folder
        reference_hvcf = args.reference_hvcf

    else:
        pangenome_hvcf = "/scratch/PHG/output/vcf_files"
        hvcf = "/scratch/PHG/output/ensambled_genomes/1740D-268-01_S1_L001_R1_001_GDB136.h.vcf"
        reference_hvcf = "/scratch/PHG/vcf_dbs/hvcf_files/MorexV3.h.vcf.gz"



    plot_folder = (f"{pangenome_hvcf}/plots/")
    hvcf_files_path = [os.path.join(pangenome_hvcf, file) for file in os.listdir(pangenome_hvcf) if file.endswith('.h.vcf.gz')]
    genomes = [os.path.basename(file).split('.')[0] for file in hvcf_files_path]
    num_genomes = len(genomes)
    hvcf_files = [os.path.basename(file) for file in hvcf_files_path]
    print(f"Your hvcf file is: \n{hvcf}\n The genomes at pangenome are currently {num_genomes}: \n{genomes}\n ")
    reference_name = os.path.basename(reference_hvcf).split('.')[0]
    print(f"Your reference hvcf file is: \n{reference_hvcf}\n The reference genome is: \n{reference_name}\n ")

    print(f"the plot folder is {plot_folder}")  

    colors_dic, ref_color = assign_colors(genomes, reference_name)

    hvcf2bed(hvcf, plot_folder, colors_dic)

    BEDSplitByChr (hvcf, plot_folder)

    sortBED (hvcf, plot_folder)

    RunPygenometracksCommandHVCF(hvcf, reference_name, plot_folder, genomes)

    AddLegendToPygenometracks (hvcf, plot_folder, genomes, reference_name)

    #Remove temporary files
    chr_list = CaptureChrFromBED(hvcf, plot_folder)
    os.remove(f"{plot_folder}TEMP_{os.path.basename(hvcf).split('.')[0]}.bed")
    os.remove(f"{plot_folder}TEMP_{os.path.basename(hvcf).split('.')[0]}_track.ini")
    os.remove(f"{plot_folder}TEMP_legend.png")
    for chr in chr_list:
        os.remove(f"{plot_folder}TEMP_{os.path.basename(hvcf).split('.')[0]}_{chr}.bed")
        os.remove(f"{plot_folder}TEMP_{os.path.basename(hvcf).split('.')[0]}_{chr}_sorted.bed")

    print("Temporary files removed")


# In[72]:


if __name__ == "__main__":
        try:
            main()    
        except KeyboardInterrupt:
            print("Script interrupted by user. Exiting...")
            sys.exit(0)
            raise KeyboardInterrupt

