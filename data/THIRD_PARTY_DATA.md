# Third-party benchmark data

The MPS files in this directory are unmodified benchmark instances from
[MIPLIB 2017](https://miplib.zib.de/). They are included to support partial
tests of the code associated with this repository. The included subset is
listed in `TEST_INSTANCES.txt` and is the union of the first 30 entries in
`prob_list_A.mat` and `prob_list_B.mat`.

## Licensing

These MPS files are third-party data and are **not** covered by the MIT
License that applies to the authors' software in this repository. Each file
remains subject to the license and attribution information recorded for its
submission by MIPLIB.

The MIPLIB 2017 paper states that submitted contributions were required to
use a license granting redistribution rights. The default license was
[Creative Commons Attribution-ShareAlike 4.0 International](https://creativecommons.org/licenses/by-sa/4.0/),
although submitters could specify a different license. Users who redistribute
or modify an instance must follow its applicable submission license and
preserve the required attribution. Authoritative instance and submission
metadata are available from the
[MIPLIB 2017 website](https://miplib.zib.de/).

The files are stored here in uncompressed MPS form for direct use by the
MATLAB/Gurobi test scripts. No changes to the mathematical models are
intended.

## Citation

Please cite the benchmark library:

> A. Gleixner, G. Hendel, G. Gamrath, et al. MIPLIB 2017: Data-Driven
> Compilation of the 6th Mixed-Integer Programming Library. *Mathematical
> Programming Computation* 13, 443-490 (2021).
> https://doi.org/10.1007/s12532-020-00194-3

The original source archives and complete instance collection are available
from the [official MIPLIB download page](https://miplib.zib.de/download).
