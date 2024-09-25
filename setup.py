from setuptools import setup, find_packages

setup(
    name='PHGv2Tools',
    version='0.1.0',
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            "CheckImputatedHaplotype = PHGv2Tools.Modules.CheckImputatedHaplotype:CheckImputatedHaplotype",
            "CheckHaplotypeAlelles = PHGv2Tools.Modules.CheckHaplotypeAlelles:CheckHaplotypeAlleles",
            "CoreRangeDetecter = PHGv2Tools.Modules.CoreRangeDetecter:PangenomeHVCFAnalizer",
            "PlotImputedHvcf = PHGv2Tools.Modules.PlotImputedHvcf:PlotHvcf",
            "PlotPangenomeChromosomesRegions = PHGv2Tools.Modules.PlotPangenomeChromosomesRegions:PlotPangenomesChromosomesRegions",
            "RangePangenomeEvolution = PHGv2Tools.Modules.RangePangenomeEvolution:RangePangenomeEvolution",
            "CondaSetup = bin.CondaSetup.sh",             
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