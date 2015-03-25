#FETCHSEQ

Extracts high quality CDS sequence from NCBI. It uses Human and mouse CCDS annotations to pull exon’s for closely related taxa. In theory, Human and Mouse CCDS coordinates are good enough to pull sequence data for any mammalian species.

Run
```
python -cds sample_input.txt -orgn Eutheria -o out_folder
```
or I would do
```
python -cds sample_input.txt -ortho Eutheria -o out_folder
```
Use of -ortho method is only recommended while extracting mammalian sequence. 

-cds takes text file with gene names [Once gene in one line]

-orgn is for selecting organismal group [Here -orgn Eutheria will pull sequences for all Eutherian mammals available in NCBI database]

-o is the output folder
