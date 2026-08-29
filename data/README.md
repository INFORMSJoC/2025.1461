# Benchmark data

The computational experiments use instances categorized under "collection"
in [MIPLIB 2017](https://miplib.zib.de/). This directory includes a compact
48-instance test subset: the union of the first 30 entries in
`prob_list_A.mat` and `prob_list_B.mat`.

The included subset is intended for `smoke_test`. Most users can test the
software immediately without downloading additional instances. Users who
want to run `test_primalFR_batch` and reproduce the complete experiments
must first download the full required dataset; the included subset is not
sufficient for that purpose.

Please cite the benchmark library as:

> A. Gleixner, G. Hendel, G. Gamrath, et al. MIPLIB 2017: Data-Driven Compilation of the 6th Mixed-Integer Programming Library. *Mathematical Programming Computation* (2021). https://doi.org/10.1007/s12532-020-00194-3

The files `prob_list_A.mat` and `prob_list_B.mat` identify all benchmark
instances and store the metadata used by the MATLAB experiment drivers.
`TEST_INSTANCES.txt` records the included test subset. These files support
partial runs of up to the first 30 cases of each experiment, for example:

```matlab
test_primalFR(2,30)
```

The complete paper experiments require additional MIPLIB files. Before
running `test_primalFR_batch`, download the full required data from the
official MIPLIB site, retain the filenames recorded in the problem lists,
and place them directly in this directory:

```text
data/
    README.md
    THIRD_PARTY_DATA.md
    TEST_INSTANCES.txt
    prob_list_A.mat
    prob_list_B.mat
    <MIPLIB instance files>
```

The MPS files are third-party data and are not covered by this repository's
MIT software license. See `THIRD_PARTY_DATA.md` for source, licensing, and
attribution information.
