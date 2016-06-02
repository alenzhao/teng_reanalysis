# Reanalysis of Teng et. al.

This is a reanalysis of the ENCODE dataset from [Teng et. al.](http://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0940-1).
The results were submitted to http://rafalab.rc.fas.harvard.edu/rnaseqbenchmark under the name `expressRF`.
In short, the only difference in the quantification pipeline was using `--rf-stranded` versus `--fr-stranded`.

NOTE: The code submitted under 'log' has an error. Please refer to the Snakefile for the correct code.
In particular, the wrong rule was copy-and-pasted.
The correct rule is `express_paper_reverse`.
