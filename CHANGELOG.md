# Changelog

## [2.1.0] - 19-05-2026
### Added
- **CLI:** `upset-plot` command to generate UpSet plots from `hapIDranges.tsv`.
- **Module:** `core-range-detector` now supports output plot formats `png`, `pdf`, and `svg` via `-f/--format`.

### Changed
- **README:** improved documentation and updated example images for plotting and analysis commands.

## [2.0.0] - 19-01-2026
### Added
- Created `CHANGELOG.md` to track version history.
- **Core CLI:** `phgtools` now uses `cli.py` with argument parsing and rich UI.
- **UI:** Added `rich` library for beautiful terminal output and help menus.
- **CLI:** Added `--version` flag to display package version.
- **CLI:** Added `--check-setup` flag to verify all dependencies are installed.
- **Module:** Added `hvcf2bed` for converting hVCF files to BED format.
- **Module:** Added `haplopainting` for visualizing haplotype blocks.
- **Module:** Added `vcf-distance` for calculating and plotting distance matrices (uses g.vcf files, outputs matrix + clustered heatmap with dendrogram).
- **Module:** Added `manage-hbed-inversions` for extracting, collapsing, summarizing, and plotting inversion blocks from `.h.bed` files.

### Changed
- Refactored all scripts to use `src` layout structure.
- Updated `pyproject.toml` to replace `setup.py`.
- **Module:** `check-imputated-haplotype` - improved data management, prints summary, calculates % in bp coverage, updated plotting.
- **Module:** `plot-imputed-hvcf` - now uses matplotlib instead of pygenometracks (faster and cleaner). Added verbose flag with summary output.
- **Module:** `check-setup` - updated to check only required dependencies (removed phg, bcftools, agc, anchorwave, tiledb, perl, pygenometracks).
- **Dependencies:** Added `seaborn` to requirements. Removed `pygenometracks` dependency.
- **README:** Complete rewrite with updated documentation for all modules.
- **EXAMPLE DATABASE:** New example database with documentation along all pipeline, creating a PHG database from scratch and to take advantage of the phgtools.

### Removed
- **Module:** `plot-pangenome-chromosomes` - functionality replaced by `haplopainting`.
- **Dependency:** Removed `pygenometracks` - no longer needed.

### Fixed
- Fixed regex matching issue in `hvcf2bed` regarding duplicate checksums.
- Fixed `ImportError` issues by standardizing the package structure.


# Previous Versions

## [1.2.0] - 2025-04-02
### Added
- **Module:** `fastaFromKey` - extract FASTA sequences from ranges using MD5 keys.
- **Module:** `CheckSetup` - validate dependencies and system requirements.
- **Module:** `RangePangenomeEvolution` - study range acquisition patterns during genome addition.
- **Module:** `GenomeIntersection` - analyze genome intersection metrics and identity percentages.
- **Module:** `CoreRangeDetector` - identify and visualize core, unique, and accessory genomic ranges.
- **Module:** `GenomeIntersectionFromMapKmers` - process PHGv2 map_kmers output for genome analysis.
- **Module:** `PlotImputedHvcf` - create ideogram visualizations of imputed h.VCF files.
- **Module:** `PlotPangenomeChromosomes` - visualize pangenome haplotypes across chromosomes.
- **Module:** `CheckHaplotypeAllelesInPangenome` - query hapIDranges.tsv files for overlapping genomic ranges.