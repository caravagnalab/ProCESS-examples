# Data generation for ENA submission

The pipeline to prepare data to be stored in ENA is the following:

1. merge `fastq.gz` for each sample and each purity level;
2. generate the manifest file for each pair of reads;
3. generate the md5check sum.

The pipeline takes in input three different parameters:

- sample identifier: in the case of tumour samples is `SPN01_1.1`, while for normal samples is the SPNID, eg `SPN01`
- sample type: either `tumour` or `normal`
- spn identifier: the SPN ID

Example of a run for a normal sample:

```
bash ena_submission.sh SPN01 normal SPN01
```

Example of a run for a tumour sample:

```
bash ena_submission.sh SPN01_1.1 tumour SPN01
```
