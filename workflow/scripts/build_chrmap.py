#!/usr/bin/env python3
# Builds a contig rename map to add/remove `chr` based on the reference .fai.
# Writes canonical mappings for 1..22,X,Y and MT/chrM. bcftools will ignore
# mappings for contigs that are not present in the VCF.

from pathlib import Path

fai_path = snakemake.input["fai"]          # provided by Snakemake
out_path = snakemake.output["map"]
log_path = snakemake.log[0]

with open(fai_path) as fh:
    first_contig = fh.readline().split("\t", 1)[0]

use_chr_prefix = first_contig.startswith("chr")

canonical = [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]
Path(Path(out_path).parent).mkdir(parents=True, exist_ok=True)

with open(out_path, "w") as fo:
    if use_chr_prefix:
        # map '1' -> 'chr1', 'MT' -> 'chrM'
        for c in canonical:
            fo.write(f"{c}\t{'chrM' if c == 'MT' else f'chr{c}'}\n")
    else:
        # map 'chr1' -> '1', 'chrM' -> 'MT'
        for c in canonical:
            fo.write(f"{'chrM' if c == 'MT' else f'chr{c}'}\t{'MT' if c == 'MT' else c}\n")

with open(log_path, "w") as lg:
    lg.write(f"Detected reference style: {'chr*' if use_chr_prefix else 'no chr prefix'} from {fai_path}\n")
    lg.write(f"Wrote mapping to {out_path}\n")
