# Info

This folder is for the variant databases. Please put the supported files in the matched folder.
You can find the links to the supportet files in the "Databases" section below.

# Usage

Run the preprocessing file manually with

```
Rscript preprocess.R
```

# Databases


GDKD: Dienstmann et al., Cancer discovery 5.2 (2015), v20.0 
==================================================================
 
==> Dowloaded from https://www.synapse.org/#!Synapse:syn2370773
==> Corresponding article: https://cancerdiscovery.aacrjournals.org/content/5/2/118.long
==> Some preprocessing done with script: "preprocess.R"
	==> Blank spaces "\s" were removed from gene names
	==> Special characters "^" and "_" were removed from columns containing PMID
	==> Variants in the same gene sharing annotations at disease, drug, evidence and association levels were aggregated into one single entry
==> Saved as GDKD.csv
   
   
CIViC: Griffith et al., Nat Genet (2017), release may differ
=====================================================================

==> Downloaded from https://civic.genome.wustl.edu/releases
==> Corresponding article: https://www.nature.com/articles/ng.3774
==> Some preprocessing done with script: "preprocess.R"
	==> Subselection of rows containing "Predictive" evidence (column "evidence_type")
	==> Variants in the same gene sharing annotations at disease, drug, evidence and association levels were aggregated into one single entry
==> Saved as CIViC.csv


TARGET: Van Allen et al., Nat Med (2014), v3  
==================================================================

==> downloaded from http://archive.broadinstitute.org/cancer/cga/target
==> Corresponding article: https://www.nature.com/articles/nm.3559
==> Manually added 24 genes from Meric-Bernstam et al., J Natl Cancer Inst.(2015)
	==> Corresponding article: https://academic.oup.com/jnci/article/107/7/djv098/913288
==> Saved as TARGET_MERIC.csv

OncoKB.csv: 
=================================================================

==> Downloaded "all actionable Variants" from http://oncokb.org/api/v1/utils/allActionableVariants.txt
==> Corresponding article: https://ascopubs.org/doi/full/10.1200/PO.17.00011
==> Preprocessing done with script "preprocess.R"
	==> Omitting rows that contain information we don't need
==> Saved as OncoKB.csv

cancer_types:
==================================================================
==> "in-house file" (from: https://github.com/jperera-bel/MTB-Report: README.md)


MTB_Databases_versions.csv: 
=================================================================
==> is produced by preprocess.R
==> Contains information about the currently used database versions/releases in CIViC.csv, GDKD.csv and TARGET_MERIC.csv
00> Is copied as additional information whenever MTB Report is used in Batch-Mode (i.e. when localmaster.r is run)
