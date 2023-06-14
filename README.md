# free_descriptors.py
A python script for generating descriptors from a directory of .mol files

### Prerequisites 
This script requires the RDKit and Mordred packages to be installed. This is best done with conda:
`conda install -c conda-forge rdkit mordred`

### Usage 
Invoke on a folder with `.mol` files to get rdkit descriptors in a CSV
__minimum example:__ 
```
python ./free_descriptors.py -i /path/to/molfiles -o output_name
```
will read all `.mol` files in `/path/to/molfiles` and create `output_name.csv` inside the same folder.

### Options
Several options exist to compute additional descriptors:
- __fragments:__ `-f` or `--fragments` flag
- __MACCS keys:__ `-M` or `--MACCS` flag
- __ECFP6 fingerprints:__ `-E` or `--ECFP6` flag
- __Mordred descriptors:__ `-m` or `--mordred` flag
- __Macrocycle descriptors:__ `-c` or `--macrocycle` flag
Optional descriptors have a tendency to take a long time to compute and to fail on some molecules. Be sure to check output files. 

### Thanks and Acknowledgements 
- This script is a modification of the rdkit_descriptors.py script by Petr Škoda 
- This script also includes code from Phyo Phyo Kyaw Zin [article](https://doi.org/10.1038/s41598-020-63192-4)
