#！/bin/bash
cd bin/
chmod 755 *
./data_norm -i ../sample_data/gene_expr.tsv -o ../sample_data/gene_expr_norm.tsv -m upqt
./data_norm -i ../sample_data/gene_expr.tsv -o ../sample_data/gene_expr_norm_hkg.tsv -m hkg -g ../sample_data/gene_list.txt
./data_filter -i ../sample_data/gene_expr_norm.tsv -o ../sample_data/gene_expr_norm_filter.tsv -p 0.75
./network_build -i ../sample_data/gene_expr.tsv -o ../sample_data/gene_co_expr.network -m spearman -p 0.001 -c 0.9 -t 2
./module_identify -i ../sample_data/gene_co_expr.network -o ../sample_data/module.txt -s 0.5 -t 2
./annotate -g ../sample_data/go-basic.obo -a ../sample_data/gene_go.assoc -n ../sample_data/gene_co_expr.network -o ../sample_data/network_go_annotation -t 2
./annotate -g ../sample_data/go-basic.obo -a ../sample_data/gene_go.assoc -m ../sample_data/module.txt -o ../sample_data/module_go_annotation -t 2
./annotate -k ../sample_data/K2ko.tsv -a ../sample_data/gene_kegg.assoc -n ../sample_data/gene_co_expr.network -o ../sample_data/network_kegg_annotation -t 2
./annotate -k ../sample_data/K2ko.tsv -a ../sample_data/gene_kegg.assoc -m ../sample_data/module.txt -o ../sample_data/module_kegg_annotation -t 2
./rwr -n ../sample_data/gene_co_expr.network -g ../sample_data/rwr_seed_genes.list -o ../sample_data/rwr_result.tsv
cd ../util/
chmod 755 *
./data_stat -i ../sample_data/gene_expr.tsv
./network_stat -i ../sample_data/gene_co_expr.network
./network_merge -i ../sample_data/test_1.network,../sample_data/test_2.network -o ../sample_data/test_merge.network
./network_extract -i ../sample_data/gene_co_expr.network -g ../sample_data/gene_list.txt -o ../sample_data/sub.network
./network_shuffle -i ../sample_data/test_1.network -o ../sample_data/random_shuffled.network
./enrich -e ../sample_data/enrichment_gene.list -b ../sample_data/background_gene.list -g ../sample_data/go-basic.obo -a ../sample_data/gene_go.assoc -p 0.05 -o ../sample_data/enrichment.go
./enrich -e ../sample_data/enrichment_gene.list -b ../sample_data/background_gene.list -k ../sample_data/K2ko.tsv -a ../sample_data/gene_kegg.assoc -p 0.05 -o ../sample_data/enrichment.kegg
./calculate_accuracy -g ../sample_data/go-basic.obo -a ../sample_data/gene_go.assoc -n ../sample_data/gene_co_expr.network -t 2
./generate_expr_matrix_from_rsem -i ../sample_data/rsem/rsem_sample.txt -o ../sample_data/rsem/rsem_gene_expr.tsv
./generate_expr_matrix_from_stringtie -i ../sample_data/stringtie/stringtie_sample.txt -o ../sample_data/stringtie/stringtie_gene_expr.tsv
./tsv_to_csv -i ../sample_data/gene_expr.tsv -o ../sample_data/gene_expr.csv
./csv_to_tsv -i ../sample_data/gene_expr.csv -o ../sample_data/gene_expr.tsv
