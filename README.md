[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc)

# A Primal Approach to Facial Reduction for SDP Relaxations of Combinatorial Optimization Problems

This archive is distributed in association with the [INFORMS Journal on Computing](https://pubsonline.informs.org/journal/ijoc) under the [MIT License](LICENSE).

The software and experiment metadata in this repository are a snapshot of the materials associated with the paper [A Primal Approach to Facial Reduction for SDP Relaxations of Combinatorial Optimization Problems](https://doi.org/10.1287/ijoc.2025.1461) by Hao Hu and Mingming Xu.

## Cite

To cite the contents of this repository, please cite both the paper and this repository, using their respective DOIs.

https://doi.org/10.1287/ijoc.2025.1461

https://doi.org/10.1287/ijoc.2025.1461.cd

Below is the BibTeX entry for citing this snapshot of the repository.

```bibtex
@misc{HuXu2026,
  author    = {Hao Hu and Mingming Xu},
  publisher = {INFORMS Journal on Computing},
  title     = {{A Primal Approach to Facial Reduction for SDP Relaxations of Combinatorial Optimization Problems}},
  year      = {2026},
  doi       = {10.1287/ijoc.2025.1461.cd},
  url       = {https://github.com/INFORMSJoC/2025.1461},
  note      = {Available for download at https://github.com/INFORMSJoC/2025.1461}
}
```

## Description

This repository contains the MATLAB implementation of the primal facial-reduction approach presented in the paper, together with the implementation of the standard facial-reduction comparison. The experiments consider DNN, Shor, and variant-Shor SDP relaxations of mixed-binary quadratic problems generated from MIPLIB 2017 instances.

The repository is a snapshot of research code used for the reported experiments rather than a general-purpose software package.

## Which test should I run?

Most users should run the partial smoke test. The repository already includes
all benchmark files needed for this test:

```matlab
addpath('scripts')
smoke_test
```

`smoke_test` runs the first five instances of experiments 2 through 5. It is
the recommended way to verify that MATLAB, YALMIP, MOSEK, Gurobi, the data,
and the experiment code work together.

The full experiment batch is intended for users reproducing the complete
computational study. Before running it, download the full required MIPLIB
dataset and place the MPS files in `data/`. The included test subset is not
sufficient for a full batch run. After installing the complete data, run:

```matlab
addpath('scripts')
test_primalFR_batch
```

## Repository layout

- `src/` contains the SDP formulations, facial-reduction procedures, and supporting MATLAB functions.
- `src/tools/` contains MATLAB linear-algebra and matrix-vectorization utilities.
- `scripts/` contains the experiment drivers and MATLAB path-configuration function.
- `data/` contains the experiment instance lists, a compact 48-instance test
  subset, and instructions for obtaining additional MIPLIB benchmark files.
- `results/` receives generated MAT result files and diary logs.

## Requirements

The experiments require:

- MATLAB;
- YALMIP;
- MOSEK;
- Gurobi; and
- the selected MIPLIB 2017 instance files.

The paper reports experiments conducted with MATLAB R2023b and Gurobi 10.0.1. The archived implementation was configured with MOSEK 11.0. MOSEK and Gurobi require separate installations and valid licenses.

If the dependencies are not already available through the local MATLAB configuration, set the following environment variables to their MATLAB package directories before running an experiment:

```matlab
setenv('YALMIP_ROOT', '/path/to/YALMIP')
setenv('MOSEK_MATLAB_ROOT', '/path/to/mosek/toolbox')
setenv('GUROBI_MATLAB_ROOT', '/path/to/gurobi/matlab')
```

The function `scripts/add_path.m` adds these directories and the repository contents to the MATLAB path. The environment-variable values are local configuration and should not be committed to the repository.

No compilation is required for the MATLAB source files.

## Data

The repository includes the union of the first 30 instances from
`data/prob_list_A.mat` and `data/prob_list_B.mat` (48 unique MPS files,
approximately 20 MB). This subset supports partial tests of every reported
experiment setting. The complete experiments require additional files from
the [MIPLIB benchmark library](https://miplib.zib.de/). See `data/README.md`,
`data/TEST_INSTANCES.txt`, and `data/THIRD_PARTY_DATA.md` for scope, source,
licensing, and attribution information.

### Partial experiment smoke test

After configuring YALMIP, MOSEK, and Gurobi, the recommended test is:

```matlab
addpath('scripts')
smoke_test
```

This runs the first five instances of experiments 2 through 5 and saves the
generated MAT files and diary logs in `results/`. The included data subset
contains every instance needed by this partial run. The instance limit does
not reduce per-instance solver limits, so experiment 5 can still take at
least 600 seconds on a selected instance.

## Replicating the experiments

The principal experiment driver is `scripts/test_primalFR.m`. Experiments 2, 3, and 4 correspond to the DNN, Shor, and variant-Shor comparisons, respectively.

Start MATLAB in the repository root and pass the desired experiment number
to the driver. Add the scripts directory once if it is not already on the
MATLAB path. For example:

```matlab
addpath('scripts')
test_primalFR(2)
```

The driver reads the corresponding problem list from `data/`, loads each
MIPLIB instance, performs the facial-reduction experiments, and saves a
MATLAB result file and diary log in `results/`.

To run only the first few instances of an experiment, pass an optional
second argument. For example, the following command runs the first three
instances of the DNN experiment:

```matlab
test_primalFR(2, 3)
```

The script `scripts/test_primalFR_batch.m` runs every configured instance of
experiments 2 through 5. It is for full reproduction and must be run only
after the complete required MIPLIB dataset has been placed in `data/`; the
included 48-instance test subset is not sufficient:

```matlab
addpath('scripts')
test_primalFR_batch
```

The implementation uses randomized computations without a fixed random seed. Consequently, independent runs can differ in individual numerical values. Runtime measurements also depend on the hardware and solver versions.

## License

The authors' software in this repository is distributed under the MIT License; see `LICENSE`. The modified third-party routine `src/tools/licols.m` is distributed under its original BSD 3-Clause license; see `src/tools/LICENSE-licols.txt` and `THIRD_PARTY_NOTICES.md`.

## Support

For assistance with the software, contact the authors using the information in `AUTHORS`.
