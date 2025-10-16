# phylostrat

Effects of stratigraphy of phylogenetic inference

## Authors

**Niklas Hohmann**\
Utrecht University\
email: n.h.hohmann [at] uu.nl\
Web page: [uu.nl/staff/NHohmann](https://www.uu.nl/staff/NHHohmann)\
Orcid: [0000-0003-1559-1838](https://orcid.org/0000-0003-1559-1838)


## Usage

Requires R, BASH, Julia, and RevBayes (https://revbayes.github.io/) as `rb`.

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

### Simulate a stratigraphic platform

In Julia's package mode, initialize the project:

```julia
Pkg> activate .
Pkg> instantiate
```
This may take some time when ran for the first time.

Then run the desired stratigraphic model, e.g.:
```julia
include("code/strat_architectures/sinusoid_run.jl")
include("code/strat_architectures/miller_2020_run.jl")
```

this will write a HDF5 file and associated `CSV`s into `data/strat/`.

To generate summary plots for both platforms in `figs/`, run:

```julia
include("code/strat_architectures/plot_section.jl")
```

## Repository structure

* code/
    * sim.R: simulate data
* data/
    * osf : Data downloaded from osf. Initially empty, filled after `code/sim.R` is run.
    * sim/: simulation data. Initially empty, filled after `code/sim.R` is run.
    * strat/: stratigraphic architectures. Initially empty, filled after the Julia code is run
* output/: results of the mcmc analyses
* logs/: console logs
* dag/: Model DAGs
* figs/: figures. Initially empty, filled with figures as the analysis is run
* LICENSE: Apache 2.0 license text
* Manifest.toml : Julia packages used
* Project.toml : Julia infrastructure
* phylostrat.Rproj : Rproject file
* snakefile_fo : run fossil only inference workflow in snakemake
* snakefile_te : run total evidence workflow in snakemake
* snakefile_sim : generate simulation data in R & store on YoDa

## License

Apache 2.0, see LICENSE file for license text.

## Funding information

Funded by the European Union (ERC, MindTheGap, StG project no 101041077). Views and opinions expressed are however those of the author(s) only and do not necessarily reflect those of the European Union or the European Research Council. Neither the European Union nor the granting authority can be held responsible for them. ![European Union and European Research Council logos](https://erc.europa.eu/sites/default/files/2023-06/LOGO_ERC-FLAG_FP.png)
