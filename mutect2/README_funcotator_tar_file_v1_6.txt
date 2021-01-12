################################################################################
# Funcotator Data Sources Package README
################################################################################

+---------------------------------------------+ 
| Data Source Version Information             |
+---------------------------------------------+ 

Version:          1.6.20190124s
Use Case:         Somatic
Source:           ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/funcotator/funcotator_dataSources.v1.6.20190124s.tar.gz
Alternate Source: gs://broad-public-datasets/funcotator/funcotator_dataSources.v1.6.20190124s.tar.gz

################################################################################

+---------------------------------------------+ 
| README                                      | 
+---------------------------------------------+ 

This is a collection of data sources to be used in conjunction with Funcotator
to annotate Germline data samples.

This folder is a top-level Data Sources Folder for The Broad Institute's 
Funcotator tool.  When running Funcotator, pass the path to this directory in
as a command-line argument:

  ./gatk Funcotator --data-sources-path PATH/TO/THIS/FOLDER ...

For more information on Funcotator, see the Funcotator tool doc or tutorial:

  https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_funcotator_Funcotator.php
  https://gatkforums.broadinstitute.org/dsde/discussion/11193/funcotator-information-and-tutorial/

For more information on GATK, see the GATK development github site:

  https://github.com/broadinstitute/gatk

################################################################################

+---------------------------------------------+ 
| Data Sources                                |
+---------------------------------------------+ 

Using this Data Sources Folder will enable the following data sources:

 achilles
--------------------
 Project Achilles is a systematic effort aimed at identifying and cataloging gene essentiality across hundreds of genomically characterized cancer cell lines.
 This data source is a mapping between gene names and cancer types.

 cancer_gene_census
--------------------
 The Cancer Gene Census (CGC) is an ongoing effort to catalogue those genes which contain mutations that have been causally implicated in cancer.

 clinvar
--------------------
 This clinvar data source is a mapping between ClinVar (a freely accessible, public archive of reports of the relationships among human variations and phenotypes, with supporting evidence) 
 and the Human Gene Mutation Database (HGMD).  This data source was published by ClinVar up until June 18, 2013.

 cosmic
--------------------
  Cosmic is the Catalogue Of Somatic Mutations In Cancer - the world's largest and most comprehensive resource for exploring the impact of somatic mutations in human cancer.

 cosmic_fusion
--------------------
  A sub-set (see cosmic) of Cosmic, containing only a mapping of Gene Names to Fusion Genes and Fusion IDs.

 cosmic_tissue
--------------------
  A sub-set (see cosmic) of Cosmic, containing only a mapping of Gene Names to mutation counts and tissue types. 

 dbsnp
--------------------
  dbSNP is world's largest database for nucleotide variations, and is part of the National Center for Biotechnology Information (NCBI), an internationally respected resource 
  for molecular biology information. As of this date, dbSNP is comprised of a large cluster of species-specific databases that contain over 12 million non-redundant sequence 
  variations (single nucleotide polymorphisms, insertion/deletions, and short tandem repeats) and over 1 billion individual genotypes from HapMap and other large-scale genotyping activities.

  This data source is a comprehensive report of short human variations formatted in VCF. It does not include genotypes, population-specific allele frequencies, or any information regarding clinical significance. It also does not include microsatellites or named variations (i.e. variations without sequence definition).

 dna_repair_genes
--------------------
  A table of genes implicated in the repair of DNA. 
  This is an update of the table cited in Wood RD, Mitchell M, & Lindahl T Mutation Research, 2005, in Science, 2001, in the reference book DNA Repair and Mutagenesis, 2nd edition, 2006, and in Nature Reviews Cancer, 2011.

 familial
--------------------
  A mapping of Familial Cancer Genes.
  This mapping goes from Gene Name to Syndrome, with Synonyms for the Syndrome itself and the Reference from which the association was created.

 gencode
--------------------
  The GENCODE project produces high quality reference gene annotation and experimental validation for human and mouse genomes.

 gencode_xhgnc
--------------------
  A mapping from GENCODE transcript ID and HGNC transcript ID.

 gencode_xrefseq
--------------------
  A mapping from GENCODE transcript ID and RefSeq RNA ID and Protein ID.

 gnomAD_exome* 
--------------------
  The Genome Aggregation Database (gnomAD), is a coalition of investigators seeking to aggregate and harmonize exome and genome sequencing data from a variety of large-scale sequencing projects, 
  and to make summary data available for the wider scientific community. In its first release, which contained exclusively exome data, it was known as the Exome Aggregation Consortium (ExAC).
  The data set provided on this website spans 125,748 exomes and 15,708 genomes from unrelated individuals sequenced as part of various disease-specific and population genetic studies, 
  totalling 141,456 individuals. Individuals known to be affected by severe pediatric disease have been removed, as well as their first-degree relatives, so this data set should serve as 
  a useful reference set of allele frequencies for severe disease studies - however, note that some individuals with severe disease may still be included in the data set, albeit likely at 
  a frequency equivalent to or lower than that seen in the general population.

	This subset of gnomAD contains only data generated from exome samples and only the allele frequency related INFO field annotations.

	By default, gnomAD_exome is disabled.  To enable it, you must extract the tar gzipped file `gnomAD_exome.tar.gz` in the data sources directory.

  * - gnomAD_exome requires a connection to the internet.

 gnomAD_genome* 
--------------------
  The Genome Aggregation Database (gnomAD), is a coalition of investigators seeking to aggregate and harmonize exome and genome sequencing data from a variety of large-scale sequencing projects, 
  and to make summary data available for the wider scientific community. In its first release, which contained exclusively exome data, it was known as the Exome Aggregation Consortium (ExAC).
  The data set provided on this website spans 125,748 exomes and 15,708 genomes from unrelated individuals sequenced as part of various disease-specific and population genetic studies, 
  totalling 141,456 individuals. Individuals known to be affected by severe pediatric disease have been removed, as well as their first-degree relatives, so this data set should serve as 
  a useful reference set of allele frequencies for severe disease studies - however, note that some individuals with severe disease may still be included in the data set, albeit likely at 
  a frequency equivalent to or lower than that seen in the general population.

	This subset of gnomAD contains only data generated from genome samples and only the allele frequency related INFO field annotations.

	By default, gnomAD_genome is disabled.  To enable it, you must extract the tar gzipped file `gnomAD_genome.tar.gz` in the data sources directory.

  * - gnomAD_genome requires a connection to the internet.

 hgnc
--------------------
  The HUGO Gene Nomenclature Committee is the only worldwide authority that assigns standardised nomenclature to human genes.  
  This is a mapping between GENCODE Gene Name and HGNC ID, along with other HGNC data (such as the HUGO gene symbol, the approved name, and more).

 oreganno
--------------------
  ORegAnno is the Open Regulatory Annotation database.
  This data source maps genome locations to ORegAnno annotations.

 simple_uniprot
--------------------
  UniProt is a comprehensive, high-quality and freely accessible resource of protein sequence and functional information.
  This is a simple mapping between Gene Name and Uniprot protein information and annotations.


