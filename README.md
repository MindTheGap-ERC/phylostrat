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

Then you can run the standard analysis using

```
bash code/inference_base.sh ID NCHAR
```

where `ID` is the replicate number (1 to 50) and `NCHAR` is the number of characters (30, 300, or 1000).

To run the analyses for scenario A and B use

```
bash code/inference_base_A.sh ID NCHAR
bash code/inference_base_B.sh ID NCHAR
```

for the skyline inference with gaps longer than 0.5 Myr in scenario A use

```
bash code/inference_skyline_A.sh ID NCHAR
```

## Repository structure

* code
    * sim.R: simulate data
* data
    * osf : Data downloaded from osf. Initially empty, filled after `code/sim.R` is run.
    * sim: simulation data. Initially empty, filled after `code/sim.R` is run.
* output: results of the mcmc analyses
