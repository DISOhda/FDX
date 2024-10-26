# FDX 2.0.1

-   Fix: kernel functions were not properly moved to internal namespace.
-   Minor documentation corrections.
-   Fixed bugs for discrete Lehmann-Romano procedure that could lead to
    incorrect results or to infinite loops.
-   Changed maintainer e-mail address.


# FDX 2.0.0

-   New features:
    -   `discrete.GR()`, `discrete.LR()`, `discrete.PB()` and their respective
        wrappers `DGR()`, `DLR()`, `DPB()`, `NDGR()` , `NDLR()` and `NDPB()` are
        now generic functions. The previously existing functionality is
        implemented in `*.default` methods.
    -   These generic functions got `*.DiscreteTestResults` methods for
        processing `DiscreteTestResults` R6 class objects from package
        `DiscreteTests` directly, so they can be used within pipes. They also
        offer some performance gains because some time-consuming checks for
        consistency are no longer necessary.
    -   For consistency of new generics and methods, the first parameter
        `raw.pvalues` needed to be renamed to `test.results`.
    -   New parameter `threshold` for `discrete.*()`, `continuous.*()`,
        `weighted.*` and their respective wrapper functions. This enables
        selection of p-values which are smaller than or equal to a certain
        value. **Note**: observed p-values and their supports are then
        re-scaled, as the p-value distributions are now becoming conditional
        distributions. If no selection is performed (i.e. `threshold = 1`),
        `print()`, `summary()` and `plot()` outputs are as before. Otherwise,
        they now respect the re-scaled conditional distributions. Additionally,
        the `FDX` S3 class output objects of these functions now include a list
        `Select` with values and information regarding selection.
    -   New parameter `pCDFlist.indices` for `discrete.*()` and their wrappers,
        which must have the same length as `pCDFlist` and may help increasing
        performance considerably. As `pCDFlist` may now include only unique
        supports, `pCDFlist.indices` must indicate the indices of the p-values
        which belong to a given support set. If `pCDFlist` has the same length
        as `test.results`, it can be omitted (by setting it to `NULL`, the
        default). If users prefer using `DiscreteTestResults` objects, they
        do not have to take care of this, as unique supports and indices are
        automatically extracted from these objects.
-   New functions `direct.discrete.GR()`, `direct.discrete.LR()` and
    `direct.discrete.PB()` as more flexible replacements for 
    `fast.discrete.GR()`, `fast.discrete.LR()` and `fast.discrete.PB()`. The
    latter have been marked as deprecated and will be removed in the future.
-   Step function evaluation in C++ code has been replaced by closely optimized
    inline functions which offer performance gains of 10-50%.


# FDX 1.0.6

-   The "fix" of 1.0.5 could cause even more trouble and has been corrected.


# FDX 1.0.5

-   Fixed a severe bug that produced completely wrong critical constants of the 
    **continuous** Guo-Romano and Lehmann-Romano procedures.


# FDX 1.0.4

-   Minor modification to `kernel_DGR_fast` C++ function to return approximated
    binomial probabilities, too.
-   Added GitHub.


# FDX 1.0.3

-   Necessary updates due to changed signatures of imported functions of package
    PoissonBinomial.


# FDX 1.0.2

-   Minor improvements to internal functions to achieve better accuracy for
    upper tail Poisson binomial probabilities.


# FDX 1.0.1

-   Fixed a bug that prevented compilation on Solaris and clang configurations.


# FDX 1.0.0

-   Initial release.
-   Added a `NEWS.md` file to track changes to the package.
