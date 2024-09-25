def PangenomeHVCFAnalizer():

    """
    This script is for analyzing the pangenome HVCF file. It will determine the number of genomes, the number of core, unique and accesion ranges and plot the results.
    """

    import os
    import matplotlib.pyplot as plt
    import argparse

    parser = argparse.ArgumentParser(description=f'{PangenomeHVCFAnalizer.__doc__}')
    parser.add_argument('-hvcf', '--hvcf-file', help='Input HVCF file', required=True)
    args = parser.parse_args()

    hvcf_file = args.hvcf_file

    if not os.path.exists(hvcf_file):
        raise ValueError('The input HVCF file does not exist')


    def NumberOfGenomes (hvcf_file):

        with open(hvcf_file, 'r') as hvcf:
            lines = hvcf.readlines()
            num_of_genomes = 0
            for line in lines:
                if line.startswith('#'):
                    continue
                else: 
                    line = line.strip()
                    line = line.split('\t')
                    keys = line[4].split(',')
                    if len(keys) > num_of_genomes:
                        num_of_genomes = len(keys)
                    #save max value of length of keys

        return num_of_genomes

    def DetermineCoreRanges(hvcf_file):

        num_of_genomes = NumberOfGenomes(hvcf_file)

        genome_numbers = list(range(1, num_of_genomes + 1))

        print(genome_numbers)

        with open(hvcf_file, 'r') as hvcf:
            lines = hvcf.readlines()
            total_ranges = 0
            core_ranges = 0
            unic_ranges = 0

            for line in lines:
                if line.startswith('#'):
                    continue
                else: 
                    total_ranges += 1
                    line = line.strip()
                    line = line.split('\t')
                    keys = line[4].split(',')
                    if len(keys) == num_of_genomes:
                        core_ranges += 1
                    if len(keys) == 1:
                        unic_ranges += 1
                        

            print("\n")
            print("here you have the whole summary:\n ")

        with open(hvcf_file, 'r') as hvcf:
            lines = hvcf.readlines()
            accesion_count = 0
            dict_ranges_count = {}
            for genome_number in genome_numbers:
                #print(f'Number of genomes: {genome_number}')
                count = 0

                for line in lines:
                    if line.startswith('#'):
                        continue
                    else:
                        line = line.split('\t')
                        keys = line[4].split(',')
                        #print(keys)
                        #print(len(keys))
                        #print(genome_number)
                        if len(keys) == genome_number:
                            count += 1

                if genome_number is not 1 and genome_number is not num_of_genomes:
                    accesion_count += count

                dict_ranges_count[genome_number] = count


                print(f'Number of ranges with {genome_number} genomes: {count}')

            print('\nTotal number of ranges: ', total_ranges)
            print('\nNumber of genomes: ', num_of_genomes)
            print('Number of core ranges: ', core_ranges)
            print("Number of accesion ranges are: ", accesion_count)
            print('Number of unique ranges: ', unic_ranges)

        #print("\nThe sum of core, unique and accesion ranges is:")
        #print(core_ranges + unic_ranges + accesion_count)
        
        if total_ranges != core_ranges + unic_ranges + accesion_count:
            raise ValueError('The sum of core, unique and accesion ranges does not match the total number of ranges')

        return total_ranges, unic_ranges, core_ranges, accesion_count, dict_ranges_count

    def plot_pangenome_ranges_intersection(hvcf_file):

        total_ranges, unic_ranges, core_ranges, accesion_count, dict_ranges_count = DetermineCoreRanges(hvcf_file)
        num_of_genomes = NumberOfGenomes(hvcf_file)

        labels = 'Core', 'Unique', 'Accesion'
        sizes = [core_ranges, unic_ranges, accesion_count]
        colors = ['gold', 'yellowgreen', 'lightcoral']
        explode = (0.1, 0, 0)  # explode 1st slice

        plt.pie(sizes, explode=explode, labels=labels, colors=colors, autopct='%1.1f%%', shadow=True, startangle=140)

        plt.axis('equal')
        plt.show()
        plt.savefig(f'{hvcf_file}_1.png')

        #get a subplot for accesion genes

        fig, ax = plt.subplots()
        ax.bar(dict_ranges_count.keys(), dict_ranges_count.values())
        ax.set_xlabel('Number of genomes')
        ax.set_ylabel('Number of ranges')
        ax.set_title('Pangenome ranges')
        #set color of first value to yellow, then green and last one red
        ax.get_children()[0].set_color('yellowgreen')
        for i in range(1, num_of_genomes):
            ax.get_children()[i].set_color('lightcoral')
        ax.get_children()[(num_of_genomes) -1].set_color('gold')
        ax.set_xticks(list(dict_ranges_count.keys()))
        ax.set_xticklabels(list(dict_ranges_count.keys()))

        plt.show()
        plt.savefig(f'{hvcf_file}_2.png')

        print(f'\nThe plots have been saved as .png files, named {hvcf_file}_1.png and {hvcf_file}_2.png')

        exit()


    plot_pangenome_ranges_intersection(hvcf_file)

PangenomeHVCFAnalizer()

