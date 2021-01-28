# Project 1 - Sequence Alignment
## Due 01/27/2021


![BuildStatus](https://github.com/MiguelGuardado/Guardado_Miguel_BMI203_HW1/workflows/HW1/badge.svg?event=push)

### Code Layout

algs.py - contains code to implement PairwiseAlignment for SmithWaterman/ NeedlemanWunsh algorithims, no main class and 
should only be used for declaring objects

main- contains all the code used for my assignment, this should be less reuseable but contains functions for sensitivity 
testing of my algorithim and can be used with correponsing files

test_align.py- unit test for my code, sort of works, I apologize for the incompletness, this message will remain as 
shame until I fix it some time in the future. 

/Images/- contain many graphs used for my assignment, can find >100 roc curves if intrested in applying differen gap_open
gap_ext pentalties on a static set of sequences.

/scoring_matrices/- contains scoring matrix that are needed for sequence alignment of amino acids

/sequences/-contains all the fasta files that were used in my assignments, can be refrences by Pospairs.txt
Negpairs.txt


### Main
Runs all code in align/\_\_main\_\_.py, useful for part 2
```
python -m align
```

### Testing
Testing is as simple as running
```
python -m pytest
```
from the root directory of this project.
