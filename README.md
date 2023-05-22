MTB_Report
=========================

This is the source code for MTB-Report. You can use this, if you want to build and start the application via command line.
NOTE: This version does not include OncoKB as part of the pipeline due to licensing issues.


MTB-Report is based on the method described in [Perera-Bel et al., 2018](https://doi.org/10.1186/s13073-018-0529-2) and uses the functions available at [MTB-Report repository](https://gitlab.gwdg.de/MedBioinf/mtb/mtb-report).

If you use this method, please cite our publication ([Perera-Bel et al., 2018](https://doi.org/10.1186/s13073-018-0529-2)).

The current release of MTB-Report uses the latest version of MTB-Report (v2.0.0), which includes the following database versions:

- [Gene Drug Knowledge Database](https://www.synapse.org/#!Synapse:syn2370773): Dienstmann et al., Cancer discovery 5.2 (2015), v20.0
- [CIViC](https://civic.genome.wustl.edu/): Griffith et al., Nat Genet (2017), release 01-Jan-2020
- [TARGET](http://archive.broadinstitute.org/cancer/cga/target): Van Allen et al., Nat Med (2014), v3
- [Meric-Bernstam et al.](https://academic.oup.com/jnci/article/107/7/djv098/913288) , J National Cancer Inst. (2015)
- [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/): Landrum et al. Nucleic Acid Res. (2018) Jan 4

The addition of ClinVar is an update included by the Bioinformatics Center at the University of Eastern Finland.

-------------------------

## Setup using Docker

MTB-Report may be run using Docker as well as directly on the host. The local deployment of MTB-Report is described in the section 'Local Setup'. Skip to this part for the Bioinformatics Centers's version of this tool.

### Download and Setup

1. To manually build the MTB-Report Docker container, clone the latest version of MTB-Report project repository: 

```
git clone https://github.com/Suhoja/MTB-Report.git
```

2. `cd` into the base directory of the project

```
cd mtb-report
```

3. Build the Docker image manually (named mtb-report):

```
docker build -t mtb-report -f ./Dockerfile .
```

  If running on a platform other than linux:

```
docker build --platform linux/x86_64/v8 -t mtb-report -f ./Dockerfile .
```

### Starting and Usage of the Application

There are two different interfaces: the web-interface and the command-line interface (CLI). This version of the tool is optimized for command-line use. Issues may arise when using the web-interface.

#### 1. The web-interface (OUT OF ORDER)

###### 1. Start the web-interface without a specific Knowledge-Database (automatically downloads and uses the CIViC-Database)
1. Start-up command:

```
docker run -d -p 3838:3838 --user shiny --name mtb-report mtb-report
```

2. Paste http://127.0.0.1:3838/mtb-report into your browser to get to the web-interface

> A user-guide for the web-interface can be found [here] (insert link when ready)


###### 2. Start the web-interface with your own Knowledge-Databases (Supports Onco-KB, CIViC, GDKD)

To include custom databases, you need to have the database directory structured in the following way: 
<pre>
databases/
|-- oncokb/
    |-- [database downloaded from oncokb]
    |-- meta.txt 
|-- civic/
    |-- [database downloaded from civic]
    |-- meta.txt 
|-- gdkd/
    |-- [database downloaded from gdkd] 
    |-- meta.txt
</pre>

Your meta.txt should contain the following lines:
<pre>
  date:    [date when the file was downloaded] 
  version: [current version of the downloaded database]
</pre>

Then run the MTB-Report Docker container with the database directory mounted in the Docker container database directory:
```
docker run -p 3838:3838 -v [database path]:/opt/MTB/databases mtb-report
```
Paste http://127.0.0.1:3838/mtb-report into your browser to get to the web-interface

-----------

#### 2. The command-line interface (RECOMMENDED)

You can run the application manually with the following command:

> MTB-Report reads files from `/work/files/`, so the data-directory `[dir]` needs to get mounted with the `-v`-flag.

```
docker run -v [dir]:/work/files/ mtb-report /work/files/[file]
```

- `[dir]`: path to the directory of your local data
- `[file]`: input-file located in `[dir]`


## Local Setup of Web Interface

### Run web application locally

To run MTB-Report on your host directly, you need to have Shiny Server installed (https://www.rstudio.com/products/shiny/download-server/). 

Configure the configuration file /etc/shiny-server/shiny-server.conf by setting the app_dir option to the MTB-Report src directory. An example shiny-server.conf is shown below. 
```
run_as shiny;

server {
  listen 3838;

  location / {

    app_dir /home/user/projects/mtb-report/src;

    log_dir /home/user/projects/mtb-report/log;

    directory_index off;
  }

```

### Run batch mode locally

To run MTB-Report directly from the command line, change into the MTB-Report base directory and run the starter script with the input file as the first script argument.
```
cd /home/user/projects/mtb-report

src/starter_c_s.sh example_data/SNV_04072022.csv
``` 

The output files are written in the same directory in the directory 'MTB_Reports'.

## Possible input types

### 1. Using a data table as input

Allowed formats: *1maf*, *vcf* and *csv*.

- <a name="maf_vcf" />For **maf** or **vcf**, only the filename is needed:

```
docker run -v [dir]:/work/files/ mtb-report /work/files/[maf | vcf]
```

- <a name="csv" />When using **csv**, additionally the mutationtype is required. (`snv`, `cnv` or `fusion`)

```
docker run -v [dir]:/work/files/ mtb-report /work/files/[csv] [snv | cnv | fusion]
```

> Example input files can be found in `example_data` (see [Testing commands with example_data](#example)).

### <a name="metadata" />2. Using a metadata file

It is also possible to use a metadata-file. This makes selecting multiple files at once possible.

```
docker run -v [dir]:/work/files/ mtb-report /work/files/[metadata-file]
```

For instructions on **writing a metadata file**:

> A simple metadata file can can be found in `example_data` (see [metadata_example.yml](#example_metadata)).


### 3. Start the command-line with your own Knowledge-Databases (Supports CIViC and GDKD. OncoKB support will be added.)

First you need a database-folder (with just the supported Knowledge-Databases you want to use) in the following structure:
<pre>
databases/ 
|-- civic/
    |-- [datbase downloaded from civic]
    |-- meta.txt 
|-- gdkd/
    |-- [datbase downloaded from gdkd] 
    |-- meta.txt
|-- oncokb/
    |-- [datbase downloaded from oncokb]
    |-- meta.txt
</pre>

Your meta.txt should contain the following lines:
<pre>
  date:    [date when the file was downloaded] 
  version: [current version of the downloaded database]
</pre>

Then run this command to start the MTB-Report Docker container.
```
docker run -v [Path databases]:/opt/MTB/databases -v [dir]:/work/files/ mtb-report /work/files/[maf | vcf]
```

## Accessing the Results

Your results are automatically added to the directory `MTB_Reports`, located in `[dir]`

There will be one metadata-file and one `csv`-file representing your results.

The new *metadata*-file consists of one `data`-part, listing your input files, and has the additional parts `tool` and `databases`

> For more detailed information about these files, go to (insert link when ready).

## <a name="example" />Testing commands with example files

This repository also has a directory `example_data` with data-examples for different formats.

Filename                                | Filetype, Mutationtype
----------------------------------------|-----
tcga_laml.maf                           | [maf-file](#maf_vcf)
M46.snps.vcf                            | [vcf-file](#maf_vcf)
Ng_et_al_2018_Supp_8_SNV_NKL.csv        | [csv, snv](#csv)
Ng_et_al_2018_Supp_2_Fusion_fh_KIJK.csv | [csv, fusion](#csv)
Ng_et_al_2018_Supp_7_CNV_My-La_1.csv    | [csv, cnv](#csv)
metadata_example.yml                    | [metadata file](#metadata)

<a name="example_metadata" />**Example metadata file:** `metadata_example.yml`

```
data:
  cancer_type: leukemia
  file_1:
    format: csv
    path: /work/files/SNV_TCGA_8.csv
    variation_type: snv
  file_2:
    format: csv
    path: /work/files/Ng_et_al_2018_Supp_7_CNV_My-La_1.csv
    variation_type: cnv
source:
  data_type: sequencing
  id: 12345
  date: 01/01/2020
```
## Changes for optimized command-line use by UEF Bioinformatics center

- Dockerfile LaTeX Linux packages modified.
- Main scripts (localmaster.r, preprocess.r) modified for .vcf input optimization.
- Helper scripts (get_druggable.r, get_levels.r, variant_annotation.r) modified for .vcf input optimization.
- Writing of the report file in .pdf format is now supported outside of the web interface. For this, the process was moved from server.r to localmaster.r and Report_local-knitr.Rnw sweave file was modified to accommodate the change.

Direct any questions regarding these modifications to Janne Suhonen (janne.suhonen@uef.fi).
