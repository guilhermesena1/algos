# My implementation of fqc


Required C++ libraries
============
 * [htslib](https://github.com/samtools/htslib) is required to process bam
   files. Please make sure you have it installed in your `lib` C++ directory.
 * [zlib](https://zlib.net) is required to read gzipped fastq files. It is
   usually installed by default in most UNIX computers and is part of the htslib
   setup, but it can also be installed with apt or brew

Installation
============
Install with the following command:
```
make all
make install
```
This will create a **bin** directory with the **fqc** executable inside, which
can either be added to your PATH variable or run locally 

Run an example as follows:
```
fqc example.fastq
```

This will generate two files : 
 * **example.fastq_qc_summary.txt** is a text file with a summary of the QC
   metrics
 * **example.fast_report.html** is the visual HTML report showing plots of the
   QC metrics summarized in the text summary. 

Copyright and License Information
=================================

Copyright (C) 2019
University of Southern California,

Current Authors: Guilherme de Sena Brandine, Andrew D. Smith

This is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your
option) any later version.

This software is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.
