# Update comments

# Tests

## Environments
* local Manjaro Linux install, R 4.1.2
* win-builder (release, oldrelease, devel)
* R-hub (configurations: "Check for CRAN", with sanitizers", "with valgrind")

## R CMD check results

### local
0 errors | 0 warnings | 1 note

Compilation used the following non-portable flag(s):
    ‘-march=x86-64’

The cause is unknown; I did not set any flags. But it does not happen on RHub
or WinBuilder.

### R-hub
0 errors | 0 warnings | 0 notes

### win-builder 
0 errors | 0 warnings | 0 notes