def PlotPangenomesChromosomesRegions():

    import os
    import gzip
    import re
    import subprocess
    import argparse

    parser = argparse.ArgumentParser(description=PlotPangenomesChromosomesRegions.__doc__)
    parser.add_argument('--hvcf_folder', "-hvcf", help='Folder with the haplotype VCF files', required=True)
    parser.add_argument('--output_folder', "-o", help='Output folder where the plot will be saved', required=True)
    parser.add_argument('--reference-hvcf', "-ref", help='Reference genome hvcf file', required=True)
    parser.add_argument('--chromosome', "-chr", help='Chromosome to plot. Enter: chrX (chr1, chr7...)', required=True)
    parser.add_argument('--region', "-r", help='Region to . Add it with a "-" dividing the start-end (100000-245000). If no value is provided, by default whole chr will be ploted', required=False)
    parser.add_argument('--reference_fasta', "-fa", help='Reference genome fasta file', required=True)

    args = parser.parse_args()

    hvcf_folder = args.hvcf_folder
    plots_folder = args.output_folder
    reference_hvcf = args.reference_hvcf
    chromosome_to_plot = args.chromosome
    region_to_plot = args.region
    ref_fasta = args.reference_fasta

    if not region_to_plot:
        region_to_plot = None

    def PreparePangenomeFilesToPlot():
        """
        Prepare the bed files and the chromosome length to plot the pangenome
        """

        def DefineHVCFs(hvcf_folder, reference_hvcf):
            """
            Define the haplotype VCF files to be plotted
            """
            haplotype_files = [f for f in os.listdir(hvcf_folder) if f.endswith('.h.vcf.gz')]
            if reference_hvcf in haplotype_files:
                haplotype_files.remove(reference_hvcf)

            return haplotype_files
        
        def DefineRegion(region_to_plot):
            """
            Define the region to plot
            """

            start = None
            end = None

            if region_to_plot:
                start, end = region_to_plot.split("-")
                output_png = f"{plots_folder}pangenome_{chromosome_to_plot}_{region_to_plot}.png"
                region_to_plot = f"{start}-{end}"

            else:
                output_png = f"{plots_folder}pangenome_FULL_{chromosome_to_plot}.png"
            
            return start, end, output_png

        def PangenomeColors():
            """
            Define the colors for the pangenome plot
            """

            haplotype_files = DefineHVCFs(hvcf_folder, reference_hvcf)
            haplotype_names = [re.sub(r'.h.vcf.gz', '', haplotype) for haplotype in haplotype_files]
            reference_name = re.sub(r'.*/|.h.vcf.gz', '', reference_hvcf)

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

            return colors_dic, ref_color, columns_header, haplotype_names, reference_name

        def CreateBedFile():
            """
            Create the bed files for each genome haplotype
            """

            colors_dic, ref_color, columns_header, haplotype_names, reference_name = PangenomeColors()

            for file in haplotype_names:
                #check if file already exists:
                print(f"starting to process {file}...")
                if os.path.exists(f"{hvcf_folder}plots/{file}_haplotype.bed"):
                    print(f"file {file}_haplotype.bed already exists, removing and making it again...")
                    os.remove(f"{hvcf_folder}plots/{file}_haplotype.bed")
                    first_line = True
                    with gzip.open(f"{hvcf_folder}{file}.h.vcf.gz", 'rt', encoding='utf-8', errors='ignore') as vcf_file, open(f"{hvcf_folder}plots/{file}_haplotype.bed", 'w') as output:
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
                    with gzip.open(f"{hvcf_folder}{file}.h.vcf.gz", 'rt', encoding='utf-8', errors='ignore') as vcf_file, open(f"{hvcf_folder}plots/{file}_haplotype.bed", 'w') as output:
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
            if os.path.exists(f"{hvcf_folder}plots/{reference_name}_haplotype.bed"):
                        print(f"file {reference_name}_haplotype.bed already exists, removing and making it again...")
                        os.remove(f"{hvcf_folder}plots/{reference_name}_haplotype.bed")
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

            return haplotype_names, reference_name

        def WholeChrLenght():

            """
            Get the whole chromosome length only if plotting the whole chromosome
            """

            reference_name = re.sub(r'.h.vcf.gz', '', reference_hvcf)

            #extract the list of real chromosomes
            with open(ref_fasta, 'r') as fasta_file:
                fasta_lines = fasta_file.readlines()
                chromosomes = [line for line in fasta_lines if line.startswith(">chr")]
                real_chromosomes = [chromosome.split(' ')[0].replace(">", "") for chromosome in chromosomes]
                print(real_chromosomes)

                chromosome_dict = {f"chr{i+1}": real_chromosomes[i] for i in range(len(real_chromosomes))}

            # Print the dictionary
            print(chromosome_dict)

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


        start, end, output_png = DefineRegion(region_to_plot)

        chromosome_length, chromosome_dict = WholeChrLenght()
       
        print("run PangenomeColors, DefineRegion, WholeChrLenght and CreateBedFile functions...")

        haplotype_names, reference_name = CreateBedFile()

        return haplotype_names, reference_name, output_png, chromosome_dict, chromosome_length

        
    
    def PlotPangenome():
        """
        Plot the pangenome
        """

        haplotype_names, reference_name, output_png, chromosome_dict, chromosome_length = PreparePangenomeFilesToPlot()

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
                track.write("line_width = 0.05\n")
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
            track.write("line_width = 0.05\n")
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

        if not region_to_plot == None:
            command_Genome_tracks = f"--tracks {track_ini} -o {output_png} --region {chromosome_dict[chromosome_to_plot]}:{region_to_plot} --dpi 300 --width 75 --fontSize 15"

        else:   #by default whole chromosome
            command_Genome_tracks = f"--tracks {track_ini} -o {output_png} --region {chromosome_dict[chromosome_to_plot]}:1-{chromosome_length} --dpi 300 --width 75 --fontSize 15 "

        print(command_Genome_tracks)

        subprocess.run(f"pyGenomeTracks {command_Genome_tracks}", shell=True)

        print(output_png)

        print(f"output file is {output_png}")
    
    PlotPangenome()

PlotPangenomesChromosomesRegions()