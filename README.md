# free_descriptors.py
A python script for generating descriptors from a directory of .mol files

### Prerequisites 
This script requires the RDKit and Mordred packages to be installed. This is best done with conda:
`conda install -c conda-forge rdkit mordred`

### Usage 
invoke on a folder with `.mol` files to get rdkit descriptors in a CSV
__minimum example:__ `python ./free_descriptors.py -i /path/to/molfiles -o output_name`
will read all `.mol` files in `/path/to/` and create `output_name.csv` inside the same folder.

### Options
several options exist to compute additional descriptors:
- fragments optional with `-f` or `--fragments` flag
- MACCS keys: -M or --MACCS flag
- ECFP6 fingerprints: -E or --ECFP6 flag
- Mordred descriptors: -m or --mordred flag
- Macrocycle descriptors: -c or --macrocycle flag
Optional descriptors have a tendency to take a long time to compute and to fail on some molecules

### Thanks and Acknowledgements 
- This script is a modification of the rdkit_descriptors.py script by Petr Škoda 
- This script also includes code from Phyo Phyo Kyaw Zin (cite)[https://doi.org/10.1038/s41598-020-63192-4]
