## Triqler

[Triqler](https://github.com/statisticalbiotechnology/triqler) is a software package that allows to compute and estimate quantitifcation error rates in label free experiments.

### Running triqler

```
$>  triqler out_triqler.tsv --decoy_pattern DECOY_ --min_samples 2 --fold_change_eval 0.5 --out proteins.tsv
```

- The file `out_triqler.tsv` is the output in the proteomicsLFQ folder for triqler
- `--decoy_pattern`: is the DECOY_ parttern selected in the proteomicsLFQ pipeline
- `--min_samples`: Min number of samples a peptide should be quantified. If this number is increased the quality of the results will be high.
- `--fold_change_eval`: log2 fold change evaluation threshold. (default: 1.0)

