from setuptools import setup, find_packages

setup(
    name='PHGv2Tools',
    version='0.1.0',
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            "check-imputated-haplotype = PHGv2Tools.Modules.CheckImputatedHaplotype:main",
            "check-haplotype-alleles = PHGv2Tools.Modules.CheckHaplotypeAllelesInPangenome:main",
            "core-range-detecter = PHGv2Tools.Modules.CoreRangeDetecter:main",
            "plot-imputed-hvcf = PHGv2Tools.Modules.PlotImputedHvcf:main",
            "plot-pangenome-chromosomes = PHGv2Tools.Modules.PlotPangenomeChromosomes:main",
            "range-pangenome-evolution = PHGv2Tools.Modules.RangePangenomeEvolution:main",
            "phgtools = PHGv2Tools.phgtools:main",
        ],
    },
    author='Joan Sarria',
    author_email='jsarria@eead.csic.es',
    description='Package to downstream analysis of pangenomes databases, working with Practical Haplotype Graph and its h.VCF files',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/jsarriaa/PHGv2Tools',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
)