# REALPHY

**REALPHY** (Reference sequence Alignment based Phylogeny builder) reconstructs phylogenetic trees of closely related bacterial (or other) genomes from short sequence reads, assembled genomes, or a mix of both. Reads and genome fragments are mapped to one or more reference genomes with Bowtie 2, polymorphic sites are extracted with quality filtering, and a tree is built with a phylogeny program of your choice.

- Website: <http://realphy.unibas.ch/realphy>
- Manual: [`doc/Realphy_manual.md`](doc/Realphy_manual.md) (PDF: <http://realphy.unibas.ch/docs/Realphy_manual.pdf>)
- Changes between versions: [`CHANGELOG.md`](CHANGELOG.md)

## Requirements

| Needed for | Program |
|---|---|
| Running REALPHY | Java 8 or newer |
| Read mapping | [Bowtie 2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) |
| Compact alignments (optional) | [samtools](http://samtools.sourceforge.net/) |
| Tree building (at least one) | [IQ-TREE](http://www.iqtree.org/), [RAxML](https://github.com/stamatak/standard-RAxML), [PhyML](http://www.atgc-montpellier.fr/phyml/binaries.php), [TREE-PUZZLE](http://www.tree-puzzle.de/) or `dnapars` from [PHYLIP](http://evolution.genetics.washington.edu/phylip/executables.html) |

The locations of the external programs are given in a `config.txt` file in the output folder; see the manual for an example.

## Building

Building from source needs a JDK 9 or newer. All dependencies are in `lib/`, so no internet access is required.

    ./build.sh

This produces a standalone `RealPhy.jar`.

## Usage

    java -jar -Xmx2000m RealPhy.jar <sequence folder> <output folder> [options]

For example, to align the sequences in `./seqs` with four threads and build the tree with IQ-TREE:

    mkdir -p out
    cp config.txt out/
    java -jar -Xmx2000m RealPhy.jar ./seqs ./out -treeBuilder 5 -thread-num 4

The sequence folder holds FASTA or GenBank files (genomes, one of which serves as reference) and FASTQ files (reads). See the manual for all options, the output folder layout and a description of the algorithm.

## Repository layout

| Path | Content |
|---|---|
| `src/` | Java sources |
| `lib/` | Bundled Java libraries (RCaller, SAM-JDK) |
| `build.sh` | Offline build script |
| `doc/` | Manual (Markdown source and PDF build script) |
| `test_run/` | Example configuration and data |

## Citation

If you use REALPHY, please cite:

> Bertels F., Silander O.K., Pachkov M., Rainey P.B., and van Nimwegen E. (2014) Automated reconstruction of whole genome phylogenies from short sequence reads. *Molecular Biology and Evolution* 31(5):1077–1088. [doi:10.1093/molbev/msu088](https://doi.org/10.1093/molbev/msu088)

GitHub shows a "Cite this repository" box generated from [`CITATION.cff`](CITATION.cff). Please also cite the programs REALPHY runs for you (Bowtie 2 and the tree builder you choose); the references are listed in the manual.

## Licence

REALPHY is free software, released under the [GNU General Public License](LICENSE), version 3 or (at your option) any later version.

The libraries in `lib/` and some source files (for example `src/util/Gamma.java`) carry their own licences; see the headers of the respective files.

## Contact

Bug reports and questions: please open an issue on GitHub.
