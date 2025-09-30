# FiLT3r

FiLT3r is a software that detects internal duplications in raw sequencing data.
FiLT3r uses an alignment-free approach, which is thus resource-frugal.

Compared to the state-of-the-art, FiLT3r is faster and uses less memory, while
having better results both for detecting the duplications and for quantifying
them.

FiLT3r was tested on capture data and RNA-Seq data for FLT3-ITD quantification. Detailed results can be found in [Baudry *et al*, 2022](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-022-04983-6). Please cite the article, if you use it.

## Compilation

```
make gatb # just once
make
```
## Docker

### Docker Hub
A docker image is available on DockerHub (`mikaels/filt3r`). You can use the `latest` tag.

### Building an own Docker image

A `Dockerfile` retrieving and compiling the program is available

You could build and run it the following way:

```shell
docker build -t filt3r:latest - < Dockerfile

# Just a filt3r -h
docker run filt3r filt3r/filt3r -h

# Example of a run
docker run -v DATA_DIR:/data -v $(pwd)/data/flt3_exon14-15.fa:/ref.fa --user $(id -u):$(id -g) filt3r filt3r/filt3r --ref /ref.fa -k 12 --sequences /data/file.fastq.gz
```

Where `DATA_DIR` is the full path to a directory storing your sequencing files and `file.fastq.gz` is the file to be analysed in that directory.
Please note that this example uses the reference for the FLT3 gene.

## Usage

`filt3r -h` provides an help on the parameters.

Note that the following parameters are **mandatory**

* `-k`: which gives the length of the *k*-mer to be used, in our paper we used
  `-k 12` with the FLT3 reference for exons 14 and 15. For other references
  you may perform some experiments to determine the best value beforehand.
* `--ref`: the reference sequence (one sequence is provided for the detection
  and quantification of FLT3/ITD in the `data/` directory).
* `--sequences` the file containing the sequences to be analysed, in case
  several files have to be analysed, they must be separated with commas
  (without spaces), *eg.* `--sequences a.fastq,b.fastq,c.fastq`.

If no output file is specified, the output filename consists of the first input filename to which a `.results.json` is appended. The output is a JSON file.

### Output file

The results are output with the JSON format.
Optionally (using `--vcf`) a VCF file can be produced.

The output file contains several fields under the JSON format :

* `nb_total`: number of reads analysed
* `nb_filtered`: number of reads matching the reference
* `nb_reads_in_events`: number of reads that overlap a breakpoint
* `percentage`: percentage of reads overlapping a breakpoint (ie. `nb_reads_in_events` / `nb_filtered`)
* `details`: provides information on each breakpoint
  * `start_pos`: Position, in the reference, at which the first *k*-mer is located after the breapoint
  * `stop_pos`: Position, in the reference, at which the last *k*-mer is located before the breapoint
  * `size`: The size of the duplication (negative when it isn't a duplication)
  * `raw_occurrence`: The number of reads where this breakpoint was detected
  * `occurrence`: The corrected number of occurrences
  * `region_coverage`: The average coverage in the region covered by the breakpoint
  * `wt_coverage`: The average coverage of the wild-type version (ie. the one where no breakpoint was detected)
  * `MT/WT`: The result of `occurrence` / `wt_coverage`
  * `VAF`: The variant allele frequency (`occurrence` / `region_coverage`)
  * `break_size`: contains the information on the size of the break of localisation. This is an internal metric that can however give clue on possible repetitions at the ends or of insertions. Break size should in theory should be equal to k-1. Deviations from this are due to the aforementioned events
  * `is_duplication`: Boolean stating whether the event is a duplication
  * `is_wt_duplication`: Boolean stating whether the event is a wild-type duplication (ie. the insertion is exactly a sequence from the reference)
  * `sequence`: the sequence of the duplication, the duplication is given in the same direction as the reference sequence. This sequence may contain some N when the duplication is not found in a single read
  * `approximate_full_sequence`: same as the `sequence` field but in that case the sequence is reconstructed at best to avoid the presence of Ns
  * `reference_pos`: The position in the full genome, if the reference sequence contains information on its coordinates in the full genome under the format `chr:start-end:strand` (eg. `13:28033760-28034429:-1`). In that case the `reference_pos` field contains the starting position in the reference of the sequence that is duplicated.


## Reproducibility

The experiments performed in the FiTL3r paper can be reproduced thanks to the
`Snakefile` under the `reproducibility/` directory. See the file
[`reproducibility/README.md`](reproducibility/README.md) for more information.
