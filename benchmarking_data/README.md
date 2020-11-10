The structure of this directory is as follows:

```
.
├── 1000_genomes
│   ├── 1000_genomes # output VCFs from 1KGBP data
│   └── 1000_genomes_pass # output VCFs from 1KBP data that had a quality value above the PASS threshold
├── dx_job_logs # logs from the DNAnexus jobs run for the benchmarking
└── hg002_benchmarks
    ├── HG002-NA24385-50x.70_percent.markdup.realigned.combined.genotyped.vcf # the output VCF from the HG002 benchmark
    └── sv_caller_outputs # individual outputs from SV callers used in the HG002 benchmark
```
