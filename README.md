# Algorithm of Hi-C data clustering and hierarchical trees comparison

The implemented algorithm, introduced by A. Rubinov (2017), builds a hierarchical tree on the given Hi-C contact matrices, and provides methods to calculate the differences between such trees.

### To run the algorithm:
Download the RubinovAlgorithm.py file, run python in the directory with the .py file and import RubinovAlgorithm class.
```
from RubinovAlgorithm import RubinovAlgorithm
```

### To read the hdf file and choose chromosome use RubinovAlgorithm class where the first argument is the path to Hi-C file, the second argument in the chromosome number:
```
ra = RubinovAlgorithm('/Users/Pushkareva/tadProject/BG3-1_20kb_IC.hdf5', 4)
```

### To visualise the tree (in the example a slice [200:400] of the Hi-C map) write:
```
ra.visualise_tree(200, 400)
```
The output is a .pdf image with the hierarchical tree visualization on the selected part of Hi-C matrix.
