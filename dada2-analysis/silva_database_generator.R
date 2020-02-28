### generates silva database for dada2

path <- "/PATH/Silva/Silva.nr_v132"
dada2:::makeTaxonomyFasta_Silva(file.path(path, "silva.nr_v132.align"), file.path(path, "silva.nr_v132.tax"), "/PATH/tax/silva_nr_v132_train_set.fa.gz")

dada2:::makeSpeciesFasta_Silva("/PATH/Silva/SILVA_132_SSURef_tax_silva.fasta.gz", "/PATH/tax/silva_species_assignment_v132.fa.gz")