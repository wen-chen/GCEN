## GCEN: an easy toolkit of Gene Co-Expression Network analysis for lncRNAs annotation  

[![platforms](https://img.shields.io/badge/platforms-Linux%20%7C%20macOS%20%7C%20Windows-yellowgreen)](https://www.biochen.org/gcen/)
[![release version](https://img.shields.io/github/v/release/wen-chen/gcen)](https://www.biochen.org/gcen)
![license GPL-3.0](https://img.shields.io/github/license/wen-chen/gcen)
![build](https://img.shields.io/travis/com/wen-chen/gcen)
[![Bioconda](https://img.shields.io/conda/dn/bioconda/gcen.svg?label=bioconda)](https://anaconda.org/bioconda/gcen)
[![star](https://img.shields.io/github/stars/wen-chen/gcen?style=social)](https://github.com/wen-chen/gcen/stargazers)
![GitHub issues](https://img.shields.io/github/issues/wen-chen/gcen)

#### Introduction  
![](https://www.biochen.org/gcen/image/workflow.png)  
GCEN is a command-line toolkit that allows biologists to easily build gene co-expression network and predict gene function, especially in RNA-Seq research or lncRNA annotation. GCEN is primarily designed to be used in lncRNA annotation, but is not limited to those scenarios. It is an efficient and easy-to-use solution that will allow everyone to perform gene co-expression network analysis without sophisticated programming skills. The recommended pipeline consists of four parts: data pretreatment, network construction, module identification, and function annotation. A README file and sample data are included in the software package. Because of its modular design, the GCEN can be easily integrated into another pipeline. Also, the multithreaded implementation of GCEN makes it fast and efficient for RNA-Seq data.

#### Highlights  
- Easy to use, without writing any code  
- High speed and low memory, can run on a personal computer  
- Suitable for RNA-Seq analysis  

#### Download  
GCEN is an open source software under the GPLv3 license. We provide source code and pre-built binaries. GCEN only supports 64-bit operating system and has been tested in Ubuntu 18.04 (Bionic Beaver), Ubuntu 20.04 (Focal Fossa), Windows 7, Windows 10, macOS 10.15 (Catalina) and macOS 11.0 (Big Sur).  

[gcen-0.5.2-linux-x86_64.tar.gz](https://www.biochen.org/gcen/download/gcen-0.5.2-linux-x86_64.tar.gz)  

[gcen-0.5.2-macOS-x86_64.tar.gz](https://www.biochen.org/gcen/download/gcen-0.5.2-macOS-x86_64.tar.gz)  

[gcen-0.5.2-windows-x86_64.zip](https://www.biochen.org/gcen/download/gcen-0.5.2-windows-x86_64.zip)  

[gcen-0.5.2-source.tar.gz](https://www.biochen.org/gcen/download/gcen-0.5.2-source.tar.gz)  

For Linux and macOS user, you can install GCEN with [conda](https://anaconda.org/bioconda/gcen). Since bioconda does not support Windows, Windows users can download the pre-built binaries directly from our website.
```
conda install gcen -c bioconda
```


#### Usage  
##### Main programs  
<details>
<summary>data_norm</summary>

```
description:
  The program data_norm normalizes the gene expression data.
usage:
  data_norm -i input_file -o output_file -m normalization_method
options:
  -i --input <input file>
  -o --output <output file>
  -m --method <upqt/median/deseq/tmm> normalization method (default: upqt)
  -v --version display GCEN version
  -h --help print help information
example:
  data_norm -i ../sample_data/gene_expr.tsv -o ../sample_data/gene_expr_norm.tsv -m upqt
```

</details>

<details>
<summary>data_filter</summary>

```
description:
  The program data_filter filter genes according to the their expression mean and standard deviation.
usage:
  data_filter -i input_file -o output_file
options:
  -i --input <input file>
  -o --output <output file>
  -c --cutoff_mean <number> mean cutoff of gene expression (default: 0.0)
  -C --cutoff_sd <number> standard deviation cutoff of gene expression (default: 0.0)
  -p --percent_mean <number> keep a proportion of total genes based mean of gene expression (default: 1.0)
  -P --percent_sd <number> keep a proportion of total genes based standard deviation of gene expression (default: 1.0)
  -v --version display GCEN version
  -h --help print help information
example:
  data_filter -i ../sample_data/gene_expr.tsv -o ../sample_data/gene_expr_filter.tsv -p 0.75
```

</details>

<details>
<summary>network_build</summary>

```
description:
  The program network_build construct gene co-expression network using gene expression matrix.
usage:
  network_build -i gene_expression_file -o co_expression_network_file
options:
  -i --input <input file>
  -o --output <output file>
  -m --method <pearson or spearman> correlation coefficient method (default: spearman)
  -l --log <log, log2 or log10> make a log(x+1) transformation (default: not transform)
  -t --thread <number> cpu cores (default: 2)
  -p --pval <number> p value cutoff (default: 0.001)
  -c --cor <number> correlation coefficient cutoff (default: 0.1)
  -s --signed <y or n> singed network (default: n)
  -f --fdr <y or n> calculate FDR (default: n)
  -v --version display GCEN version
  -h --help print help information
example:
  network_build -i ../sample_data/gene_expr.tsv -o ../sample_data/gene_co_expr.network -m spearman -p 0.001 -f y
```

</details>

<details>
<summary>module_identify</summary>

```
description:
  The program identify the gene modules using the gene co-expression network.
usage:
  module_identify -i input_file -o output_file
options:
  -i --input <input file>
  -o --output <output file>
  -s --similarity <number> similarity cutoff (default: 0.5)
  -t --thread <number> cpu cores (default: 2)
  -v --version display GCEN version
  -h --help print help information
example:
  module_identify -i ../sample_data/gene_co_expr.network -o ../sample_data/module.txt
```

</details>

<details>
<summary>annotate</summary>

```
description:
  The program annotate can perform GO/KEGG annotation based on network or module.
usage:
  annotate -g go-basic.obo -a gene_go_association_file -n input_network -o out_dir
options:
  -g --go <go-basic.obo file>
  -k --kegg <kegg information> (if the -g/--go is specified, the -k/--kegg are ignored)
  -a --assoc <gene/go association file>
  -n --network <network file>
  -m --module <module file> (if -n is specified, the -m is ignored)
  -p --pval <number> p value cutoff (default: 0.05)
  -o --output <output directory>
  -t --thread <number> cpu cores (default: 2)
  -v --version display GCEN version
  -h --help print help information
examples:
  ./annotate -g ../sample_data/go-basic.obo -a ../sample_data/gene_go.assoc -n ../sample_data/gene_co_expr.network -o ../sample_data/network_go_annotation
  ./annotate -g ../sample_data/go-basic.obo -a ../sample_data/gene_go.assoc -m ../sample_data/module.txt -o ../sample_data/module_go_annotation
  ./annotate -k ../sample_data/K2ko.tsv -a ../sample_data/gene_kegg.assoc -n ../sample_data/gene_co_expr.network -o ../sample_data/network_kegg_annotation
  ./annotate -k ../sample_data/K2ko.tsv -a ../sample_data/gene_kegg.assoc -m ./sample_data/module.txt -o ../sample_data/module_kegg_annotation
```

</details>

<details>
<summary>rwr</summary>

```
description:
  The program rwr can predict potential funcation associated genes based on RWR (Random Walk with Restart) algorithm.
usage:
  rwr -n input_network -g gene_list -o output_result
options:
  -n --network <network file>
  -g --gene <seed genes list file>
  -r --gamma <number> restart probability (default: 0.5)
  -p --prob <number> probability cutoff (defalut: 0.01)
  -o --output <output file>
  -d --directed_network the input network is directed (defalut: the input network is undirected)
  -w --weighted_network the edge weights of network will be considered (defalut: all edge weights of network are set to 1.0)
  -W --weighted_gene the weights of seed genes will be considered (defalut: all weights of seed genes are set to 1.0)
  -v --version display GCEN version
  -h --help print help information
example:
  rwr -n ../sample_data/gene_co_expr.network -g ../sample_data/rwr_seed_genes.list -o ../sample_data/rwr_ranked_gene.tsv
```

</details>

##### Utilities  

<details>
<summary>data_stat</summary>

```
description:
  The program data_stat calculate the statistics of gene expression matrix.
usage:
  data_stat -i input_file
options:
  -i --input <input file>
  -v --version display GCEN version
  -h --help print help information
example:
  data_stat -i ../sample_data/gene_expr.tsv
```

</details>

<details>
<summary>network_merge</summary>

```
description:
  The program network_merge merge two or more networks.
usage:
  network_merge -i input_files -o output_file
options:
  -i --input <input files> multiple files are separated by commas
  -o --output <output file>
  -c --cor <number> correlation coefficient cutoff (default: 0.5)
  -h --help print help information
example:
  network_merge -i ../sample_data/test_1.network,../sample_data/test_2.network -o ../sample_data/test_merge.network
```

</details>

<details>
<summary>enrich</summary>

```
description:
  The program enrich can perform GO/KEGG enrichment.
usage:
  enrich -e enrichment_gene_list_file -b background_gene_list_file -g go-basic.obo -a gene_go_association_file -p p_value_cutoff -o out_put_file
options:
  -e --enrich <enrichment gene list file>
  -b --background <background gene list file>
  -g --go <go-basic.obo file>
  -k --kegg <kegg information> (if the -g/--go is specified, the -k/--kegg are ignored)
  -a --assoc <gene/go association file>
  -p --pval <number> p value cutoff (default: 0.05)
  -o --output <output file>
  -v --version display GCEN version
  -h --help print help information
examples:
  enrich -e ../sample_data/enrichment_gene.list -b ../sample_data/background_gene.list -g ../sample_data/go-basic.obo -a ../sample_data/gene_go.assoc -p 0.05 -o ../sample_data/enrichment.go
  enrich -e ../sample_data/enrichment_gene.list -b ../sample_data/background_gene.list -k ../sample_data/K2ko.tsv -a ../sample_data/gene_kegg.assoc -p 0.05 -o ../sample_data/enrichment.kegg
```

</details>

<details>
<summary>generate_expr_matrix_from_rsem</summary>

```
description:
  The program generate_expr_matrix_from_rsem generate gene expression matrix from RSEM outputs.
usage:
  generate_expr_matrix_from_rsem -i input_file -o output_file
options:
  -i --input <input file> a text file with sample ID and path to its RSEM result file on each line
  -o --output <output file>
  -t --tpm output TPM value instead of FPKM vaule
  -v --version display GCEN version
  -h --help print help information
example:
  generate_expr_matrix_from_rsem -i ../sample_data/rsem/rsem_sample.txt -o ../sample_data/rsem/rsem_gene_expr.tsv
```

</details>

<details>
<summary>generate_expr_matrix_from_stringtie</summary>

```
description:
  The program generate_expr_matrix_from_stringtie generate gene expression matrix from StringTie outputs.
usage:
  generate_expr_matrix_from_stringtie -i input_file -o output_file
options:
  -i --input <input file> a text file with sample ID and path to its GTF file on each line
  -o --output <output file>
  -t --tpm output TMP value instead of FPKM vaule
  -v --version display GCEN version
  -h --help print help information
example:
  generate_expr_matrix_from_stringtie -i ../sample_data/stringtie/stringtie_sample.txt -o ../sample_data/stringtie/stringtie_gene_expr.tsv
```

</details>


#### Recommended pipeline  

<details>
<summary>Step 1: data pretreatment
</summary>

```
./data_norm -i ../sample_data/gene_expr.tsv -o ../sample_data/gene_expr_norm.tsv -m upqt
./data_filter -i ../sample_data/gene_expr_norm.tsv -o ../sample_data/gene_expr_norm_filter.tsv -p 0.75
```

</details>

<details>
<summary>Step 2: co-expression network construction
</summary>

```
./network_build -i ../sample_data/gene_expr_norm_filter.tsv -o ../sample_data/gene_co_expr.network -m spearman -p 0.001 -c 0.8 -f y -t 6
```

</details>

<details>
<summary>Step 3: module identification (optional)</summary>

```
./module_identify -i ../sample_data/gene_co_expr.network -o ../sample_data/module.txt -s 0.5 -t 6
```

</details>

<details>
<summary>Step 4: function annotation</summary>

```
# network based annotation
./annotate -g ../sample_data/go-basic.obo -a ../sample_data/gene_go.assoc -n ../sample_data/gene_co_expr.network -o ../sample_data/network_go_annotation
./annotate -k ../sample_data/K2ko.tsv -a ../sample_data/gene_kegg.assoc -n ../sample_data/gene_co_expr.network -o ../sample_data/network_kegg_annotation

# module based annotation (optional)
./annotate -g ../sample_data/go-basic.obo -a ../sample_data/gene_go.assoc -m ../sample_data/module.txt -o ../sample_data/module_go_annotation
./annotate -k ../sample_data/K2ko.tsv -a ../sample_data/gene_kegg.assoc -m ../sample_data/module.txt -o ../sample_data/module_kegg_annotation

# identify genes with specific functions based on RWR (optional)
rwr -n ../sample_data/gene_co_expr.network -g ../sample_data/rwr_interested_gene.list -o ../sample_data/rwr_result.tsv
```

</details>


#### Data format  
To understand the format of the input and output files for each program, please take a look at the sample data included in the software package.  

#### Contact  
We appreciate your interest in our work. For further information or if you have any questions, please do not hesitate to contact us (Wen Harold Chen, chenwen@biochen.org).
