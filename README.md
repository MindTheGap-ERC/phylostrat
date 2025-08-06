# phylostrat

Effects of stratigraphy of phylogenetic inference

## Usage

Requires R, BASH, and RevBayes (https://revbayes.github.io/) as `rb`.

### Generate simulation data

To simulate the data, run

```
Rscript code/sim.R
```

in BASH. This will download approx. 190 Mb of data from OSF if needed, and then generate all simulation data under `data/sim/`. This will run a few minutes and produce approx 50 Mb of data.

The repository contains multiple types of analyses:
* fossil only, suffix `fo`
* total evidence, suffix `te`
* constant rate FBD analyses, suffix `const` or none
* skyline FBD inference with sampling during destructive intervals set to 0, suffix `skyline_A`
* skyline FBD inference with all sampling parameters estimates, suffix `skyline_A_full`

## Repository structure

* code
    * sim.R: simulate data
* data
    * osf : Data downloaded from osf. Initially empty, filled after `code/sim.R` is run.
    * sim: simulation data. Initially empty, filled after `code/sim.R` is run.
* output: results of the mcmc analyses
* logs: console logs
* dag: Model DAGs
* LICENSE: Apache 2.0 license text
* phylostrat.Rproj : Rproject file
* snakefile_fo : run fossil only inference workflow in snakemake
* snakefile_te : run total evidence workflow in snakemake
* snakefile_sim : generate simulation data in R & store on YoDa