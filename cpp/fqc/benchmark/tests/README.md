Placeholder directory for test datasets
========================================
This directory can be used to download the SRR files. With the SRA toolkit
installed, users can run the following command to download the benchmark
datasets:

```
fastq-dump --outdir . --gzip --skip-technical --readids --read-filter pass
--dumpbase --split-3 --clip SRRXXXXXX
```

Where **SRRXXXXXX** should be replaced by each of the following SRRs used in the
paper:

 - SRR10124060
 - SRR10143153
 - SRR3897196
 - SRR9624732
 - SRR1853178
 - SRR6387347
 - SRR1772703
 - SRR9878537
 - SRR6059706

The human genome nanopore file can be downloaded from the [Human Whole Genome
Sequencing
Project](https://github.com/nanopore-wgs-consortium/NA12878/blob/master/nanopore-human-genome/rel_3_4.md) (file [FAB49164](http://s3.amazonaws.com/nanopore-human-wgs/rel3-nanopore-wgs-4045668814-FAB49164.fastq.gz))


