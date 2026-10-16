# Changelog

## 1.14

### Added
- Parallel execution of the alignment step. The new option `-thread-num [integer]` (default 1) sets the number of sequence files that are fragmented and aligned by bowtie2 at the same time. The number is also passed to IQ-TREE.
- IQ-TREE as tree building method (`-treeBuilder 5`). The executable is configured with the key `IQTree` in `config.txt` (default `iqtree`).
- Option `-alignment-only` to stop after the alignment step.
- Command line options for IQ-TREE (and the other tree programs) can be read from the file given with `-treeOptions`.
- `build.sh` builds a standalone `RealPhy.jar` offline from the sources and the jars in `lib/`.
- Manual in Markdown in `doc/` with a script to build the PDF.

### Changed
- PhyML output file `*_phyml_tree` (new PhyML versions) is recognised in addition to `*_phyml_tree.txt`.
- Java 8 or newer is required (previously the manual stated Java 1.6).
- Version number is 1.14 (was 1.13).
- Manual and `-help` text corrected to match the program (defaults of `-ref`, `-root` and `-config`, `-refN`, `-clean`, `-merge`, `-gaps`, `-gapThreshold`, `-varOnly`, `-d`; option names `-polyThreshold` and `-treeBuilder`) and `-ref-list` documented.
