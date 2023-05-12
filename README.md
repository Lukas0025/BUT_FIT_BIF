# BioPython for create ancestrals seqvences
#### Semestral project to BIF on BUT FIT

Create ancestrals amino acides seqvences for Phylogenetic tree

## Used packages

* BioPython
* Pandas

## Run

just start `main.py` using python3 interpreth

```sh
python3 main.py
```

## Need files

for run program prgram need some files. This files can be chaged in top of python script.

```python
TREE_FILE       = "tree.tre"
MSA_FILE        = "msa.fasta"
ANCESTRALS_FILE = "ancestrals.csv"
OUT_DIR         = "output"
```

* `TREE_FILE` contains fylogenetic tree in Newick (with confidence)
* `MSA_FILE` cantains aligmented leaves seqvences in FASTA 
* `ANCESTRALS_FILE` contains posterior probabilty of every aminoacid in every seqvence