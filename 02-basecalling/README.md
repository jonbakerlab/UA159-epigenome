# 02. Basecalling `.pod5` files using `dorado basecaller`
`dorado basecaller` basecalls the POD5 files, which contains the raw electrical signals measured by the nanopore sequencer, and outputs it into a **single** `.bam` file format. Note that **`dorado basecaller` utilizes GPUs to run.** With the `--modified-bases-models` flag, the `.bam` file will contain MM/ML values of modified bases of the provided models known as a `modBAM` file.

Assuming multiple samples were run on single flow cell run, all samples are multiplexed in the single `.bam` output file and therefore must be demultiplexed using `dorado demux`. Here, the sequencing library prep kit name (SQK-NBD114-24) must be provided for the `--kit-name` argument to identify the unique barcodes and demultiplex the samples based on that. The basecalling pipeline can be found in the `sbatch_nativeDNA_basecall_mods.sh` script.

## Additional information
### Troubleshooting `dorado`
`dorado` can struggle with memory allocation when running the program. Meaning that if too much memory is provided, it can over utilize the memory quickly and run into an error. Through trial and error, the parameters above seem to a good standard for minimizing total run-time while not running into an error. However, if the total size of the `pod5` files supplied is very large, the general rule of thumb would be to either decrease the total memory available or lowering the batch size via the `--batch-size` flag. In addition, if the run does run into an error, the flag `--resume-from` can be used with the argument being the path to the incomplete modBAM file. Note that if this is used, the output file **must** have a different name than the argument provided for `--resume-from`, e.g. `multiplexed-2.bam`. This can be implimented as shown below:

```sh
path/to/dorado-1.4.0-linux-x64/bin/dorado basecaller \
path/to/dorado-1.4.0-linux-x64/bin/models/dna_r10.4.1_e8.2_400bps_sup@v5.2.0 \
path/to/pod5_directory \
--modified-bases-models path/to/dorado-1.4.0-linux-x64/bin/models/dna_r10.4.1_e8.2_400bps_sup@v5.2.0_6mA@v1,path/to/dorado-1.4.0-linux-x64/bin/models/dna_r10.4.1_e8.2_400bps_sup@v5.2.0_4mC_5mC@v1 \
--modified-bases-batchsize 128 \
--min-qscore 9 \
--batchsize 128 \
--resume-from ./dorado-output/multiplexed.bam \
> ./dorado-output/multiplexed-2.bam
```
### Outsourcing nanopore sequencing to Plasmidsaurus
Plasmidsaurus is a whole genome sequencing company that sequences samples using Oxford Nanopore Technologies. Typically they only only provide `.fastq` files from their services, however upon request, they will provide a raw `.pod5` file, which is required to basecall and identify the modified bases. The `sbatch_DNA_modBAM_plasmidsaurus.sh` scripts can accomodate `.pod5` directly sourced from Plasmidsaurus.