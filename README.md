# ABraOM VEP Plugin

This Ensembl VEP plugin adds transcript-level population frequency annotations from the ABraOM database.

## Overview

The plugin queries a tabix-indexed ABraOM file and adds selected database columns as VEP annotations. By default, it exports only the `Frequencies` column as `abraom_freq`, but you can specify additional columns to include.

## Installation

1. Copy `ABraOM.pm` to the VEP plugins directory, for example:

```bash
mv ABraOM.pm ~/.vep/Plugins
```

2. Download the ABraOM data file and prepare it for tabix:

```bash
tar -xvzf SABE1171.Abraom.clean.tar.gz
cat SABE1171.Abraom.clean.tar.gz | (head -n 1 && tail -n +2 | sort -t$'\t' -k 1,1 -k 2,2n) > SABE1171.Abraom_sorted.tsv
sed '1s/.*/#&/' SABE1171.Abraom_sorted.tsv > SABE1171.Abraom_final.tsv
bgzip SABE1171.Abraom_final.tsv
tabix -f -s 1 -b 2 -e 2 SABE1171.Abraom_final.tsv.gz
```

## Usage

### Annotate only the default frequency

```bash
vep -i variations.vcf --plugin ABraOM,file=/path/to/abraom/SABE1171.Abraom_final.tsv.gz
```

This adds the `abraom_freq` annotation to the output.

### Annotate specific columns

```bash
vep -i variations.vcf --plugin ABraOM,file=/path/to/abraom/SABE1171.Abraom_final.tsv.gz,cols=Frequencies,Gene.refGene,PredictedFunc.refGene
```

This adds columns with the `abraom_` prefix, such as `abraom_Gene.refGene`.

## Plugin parameters

- `file=`: path to the tabix-indexed ABraOM file (`.tsv.gz`). Required.
- `cols=`: comma-separated list of columns to export. Default: `Frequencies`.

## Supported columns

The plugin maps ABraOM file columns to VEP annotations. Supported column examples include:

- `Frequencies`
- `Gene.refGene`
- `PredictedFunc.refGene`
- `avsnp150`
- `FILTER`
- `CEGH_Filter`
- `HomozygousALT_count`
- `Hemizygous_count`
- `Allele_number`
- `Allele_ALT_count`
- `Cohort`

## Requirements

- Ensembl VEP with plugin support
- Perl
- `tabix` and `bgzip`
- ABraOM database in tabix-indexed format

## License

Distributed under the `Apache License 2.0`. See the header of `ABraOM.pm` for more details.

## Contact

- Email: romualdofarm@gmail.com
