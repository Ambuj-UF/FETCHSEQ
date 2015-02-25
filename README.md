#FETCHSEQ

Extracts high quality CDS sequence from NCBI. It uses Human and mouse CCDS annotations to pull exon’s for closely related taxa. In theory, Human and Mouse CCDS coordinates are good enough to pull sequence data for any mammalian species.

Run
```
python -cds sample_input.txt -orgn Eutheria -o out_folder
```