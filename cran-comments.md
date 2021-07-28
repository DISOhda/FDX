# Update comments

# Tests

## Environments
* local Windows 10 Education (20H2) install, R 4.0.3
* win-builder (release, oldrelease, devel)
* R-hub (configurations: "with sanitizers", "with valgrind")

## R CMD check results

### local
0 errors | 0 warnings | 1 note

* checking compiled code ... NOTE
  Note: information on .o files for i386 is not available
  Note: information on .o files for x64 is not available
  File 'D:/Nextcloud/R/PoissonBinomial.Rcheck/PoissonBinomial/libs/i386/PoissonBinomial.dll':
    Found 'abort', possibly from 'abort' (C), 'runtime' (Fortran)
    Found 'exit', possibly from 'exit' (C), 'stop' (Fortran)
    Found 'printf', possibly from 'printf' (C)
  File 'D:/Nextcloud/R/PoissonBinomial.Rcheck/PoissonBinomial/libs/x64/PoissonBinomial.dll':
    Found 'abort', possibly from 'abort' (C), 'runtime' (Fortran)
    Found 'exit', possibly from 'exit' (C), 'stop' (Fortran)
    Found 'printf', possibly from 'printf' (C)

The cause is unknown; none of the mentioned functions is used. This NOTE did not
appear in any of the previous versions or in any online check (see below).

### R-hub
0 errors | 0 warnings | 0 notes

### win-builder 
0 errors | 0 warnings | 0 notes