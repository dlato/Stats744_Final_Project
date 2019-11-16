There are two papers that are related to the data:
1) substitutions_paper_draft.pdf
2) gene_expression_paper_draft.pdf

The substitutions paper has information about the substitutions data (or mutations) and the selection analysis (if a gene is under positive or negative selection)
The substitutions paper has some of the statistical methods that are out of date, but the results in the supplemental file, figures, and tables are all correct and up to date.

The gene expression paper is complete and up to date, including the supplementary file.

The Data files:
1) Substitutions data files: (*_bidirectional_data.csv)
  - each file is set up the same way and there is one file per bacterial replicon:
    - ecoli = Escherichia coli
    - bass = Bacillus subtilis
    - strep = Streptomyces
    - sinoC = Sinorhizobium meliloti Chromosome
    - pSymA = Sinorhizobium meliloti pSymA (2nd chromosome)
    - pSymB = Sinorhizobium meliloti pSymB (3rd chromosome)
  - the columns are as follows from left to right:
    1) row number (leftover from me saving the data frame from R) (useless)
    2) "branch" : metadata, (related to the phylogenetic analysis, complicated but irrelevant for what we want to do)
    3) "from1" : metadata, node 1 from branch (again related to the phylogenetic cnalysis, irrelevant for what we want to do)
    4) "to1" : metadata, node 2 from branch (again related to the phylogenetic cnalysis, irrelevant for what we want to do)
    5) "position" : this is the ORIGINAL genomic position, we do NOT want to use this one!
    6) "change" : is there a substitution at this site? 1 = yes, 0 = no
    7) "from2" : metadata, if there was a substitution, it changed from this column to the next coloumn ("to2"). for example no substitution would be NA, a substitution would be "A", "G", "C", or "T")
    8) "to2" : metadata, if there was a substitution, it changed from the previous column ("from2") to this column. for example no substitution would be NA, a substitution would be "A", "G", "C", or "T")
    9) "prob" : metadata, probability of the change/substitution
    10) "bidir_pos" : position from origin of replication, WE WANT TO USE THIS ONE as the genomic position!!


2) Gene Expression: (data/*_gene_exp_bidirectional_data.csv)
  - each file is set up the same way and there is one file per bacterial replicon:
    - ecoli = Escherichia coli
    - bass = Bacillus subtilis
    - strep = Streptomyces
    - sinoC = Sinorhizobium meliloti Chromosome
    - pSymA = Sinorhizobium meliloti pSymA (2nd chromosome)
    - pSymB = Sinorhizobium meliloti pSymB (3rd chromosome)
  - the columns are as follows from left to right:
    1) row number (leftover from me saving the data frame from R) (useless)
    2) "X" : row number (leftover from me saving the data frame from R) (useless)
    3) "Id" : gene id for each gene in the genome
    4) "Start" : metadata, ORIGINAL genomic start position of this gene
    5) "End" : metadata, ORIGINAL genomic end position of this gene
    6) "tmp_pos" : position in the genome from the origin of replication, WE WANT TO USE THIS ONE AS THE GENOMIC POSITION!!
    7) "Exp" : expression value for that gene in normalized CPM
    
3) Selection Data: (data/*_selection_data.csv)
  - each file is set up the same way and there is one file per bacterial replicon:
    - ecoli = Escherichia coli
    - bass = Bacillus subtilis
    - strep = Streptomyces
    - sinoC = Sinorhizobium meliloti Chromosome
    - pSymA = Sinorhizobium meliloti pSymA (2nd chromosome)
    - pSymB = Sinorhizobium meliloti pSymB (3rd chromosome)
  - the columns are as follows from left to right:
    1) "gene_name" : Gene Id
    2) "dN" : nonsynonymous mutation rate, i.e. mutations that DO change the amino acid and SHOULD alter the function of the gene
    3) "dS" : synonymous mutation rate, i.e. mutations that do not change the amino acid and should not be changing the functionality of the gene
    4) "omega" : ratio of nonsynonymous / synonymous mutation rates, omega > 1 means that there is positive selection on that gene, omega < 1 mean there is negative selection on that gene, omega = 1 (or close to it) means there is neutral selection on the gene
    5) "midpoint" : metadata, ORIGINAL genomic end position of this gene
    6) "tmp_pos" : position in the genome from the origin of replication, WE WANT TO USE THIS ONE AS THE GENOMIC POSITION!!
    
    
    
