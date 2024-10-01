import argparse
from PHGv2Tools.Modules.CheckHaplotypeAllelesInPangenome import main as check_haplotype_alleles_main
from PHGv2Tools.Modules.CheckImputatedHaplotype import main as check_imputated_haplotype_main
from PHGv2Tools.Modules.CoreRangeDetecter import main as core_range_detecter_main
from PHGv2Tools.Modules.PlotImputedHvcf import main as plot_imputed_hvcf_main
from PHGv2Tools.Modules.PlotPangenomeChromosomes import main as plot_pangenome_chromosomes_main
from PHGv2Tools.Modules.RangePangenomeEvolution import main as range_pangenome_evolution_main

def main():
    parser = argparse.ArgumentParser(prog='phgtools', description='PHG Tools Command Line Interface')
    subparsers = parser.add_subparsers(dest='command')

    # Add subcommands
    subparsers.add_parser('check-haplotype-alleles', help='Check haplotype alleles in pangenome')
    subparsers.add_parser('check-imputated-haplotype', help='Check imputated haplotype')
    subparsers.add_parser('core-range-detecter', help='Core range detecter')
    subparsers.add_parser('plot-imputed-hvcf', help='Plot imputed HVCF')
    subparsers.add_parser('plot-pangenome-chromosomes', help='Plot pangenome chromosomes')
    subparsers.add_parser('range-pangenome-evolution', help='Range pangenome evolution')

    args = parser.parse_args()

    if args.command == 'check-haplotype-alleles':
        check_haplotype_alleles_main()
    elif args.command == 'check-imputated-haplotype':
        check_imputated_haplotype_main()
    elif args.command == 'core-range-detecter':
        core_range_detecter_main()
    elif args.command == 'plot-imputed-hvcf':
        plot_imputed_hvcf_main()
    elif args.command == 'plot-pangenome-chromosomes':
        plot_pangenome_chromosomes_main()
    elif args.command == 'range-pangenome-evolution':
        range_pangenome_evolution_main()
    else:
        parser.print_help()

if __name__ == '__main__':
    main()