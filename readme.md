# DMS Missense Variants v2

Update on designing missense variant primers across windows of interest

**STATUS:** 
- use bool logic for synonymous variants and remove stop codons in mut windows
- need to check the final codon table
- need a relative path setup for the dms_codon_table_v2

## Input:
- Wildtype .gb file
    - accomodate reverse complement?
- Options:
    - Vector .gb file
    - synonymous mutation rate (def .05)
    - random seed (def 42)
    - codon table (def 'Standard')
    - homology arm length (def 20)
    - primer anneal temp (def 50)
    - max primer length (aim for 60bp)
    
## Output:
- .tsv codon table
    - primer name
    - primer seq
    - codon variant (e.g. ACC71ACB)
- fasta option?

## Updates:
- missense table correction
- more complete synonymous variants table (random seed choice)
- random seed for synonymous variants as controls (reproducible)
- .gb input file
    - reverse complement
    - provide window names

## RESOURCES
- [Organizing python package imports](https://towardsdatascience.com/whats-init-for-me-d70a312da583)

## Sketch

imports
arguments

generate missense codon dictionary (script) # UPDATE: use selected dms v2 codon table
generate synonymous codon dictionary (script) # UPDATE: use selected dms v2 codon table

parse .gb file, loop for each matching feature "window"
    - check that window is divisible by 3, codons
    - check for upstream homology and downstream primer space (20bp, 40bp)

    assign sub-window start index value
    define which codons will contain synonymous controls at 5% frequency (based on wt)
    define synonymous codons (based on vector)
    (create all the codon variants immediately, then define subwindows based on primers)
    
    begin sub-window loop:
        homology arm
        primer design
            - redefine sub-window start index
        drop-in missense sub-window
        
aggregate into dataframe