# CRISPRutils
Utility scripts for tedious everyday tasks

## Installation
```
# Download CRISPRutils
git clone https://github.com/CRISPRlab/CRISPRutils.git

# Modify PATH variable (Mac, Linux) | (.bash_profile, .profile, .bashrc, etc.)
vi ~/.bash_profile

# Adjust PATH variable, then save
export PATH=$PATH:/path/to/cloned/repo/CRISPRutils

# Load newly modified PATH variable
source ~/.bash_profile
```

---
# Scripts
### < append_spacers.sh >
A simple script that locates all files with the extension '*_spacers.fa*' in the **current directory**, then appends them into a single FASTA file, named *spacers_concat.fasta*.

#### Software Requirements
- **Bash** >= 3.2.57

#### Usage:

```
cd directory_with_spacer_files
append_spacers.sh
```


---
### < blast_parser.py >
Parses returned BLAST hits into CSV format, and assists in PAM prediction by fetching flanking regions from [Entrez](https://www.ncbi.nlm.nih.gov/Class/MLACourse/Original8Hour/Entrez/) by each BLAST hit's Accession number.

#### Software Requirements
- **Python** >= 2.7 or 3 AND **pip** (install pip for Python 2 and pip3 for Python 3)
  - Mac - using homebrew:

   `brew install python`

   `sudo easy_install pip`
  - Ubuntu:

   `sudo apt-get install -y python3-pip`

- **Biopython**

  `pip install biopython`
  or

  `pip3 install biopython`

- **Pandas**

  `pip install pandas`

**Important Note:** when using this script to predict the PAM sequence (-p option), we use Entrez to access genomes via Accession number. Please read [NCBI's Entrez user requirements](https://www.ncbi.nlm.nih.gov/books/NBK25497/#chapter2.Usage_Guidelines_and_Requiremen) prior to use. At minimum, please modify **blast_parser.py at line 155** with the addition of your **email address**. NCBI will attempt to contact you before killing (or banning) your jobs in cases of overuse or abuse. Otherwise, your job(s) could be killed without warning. Providing a fake email address is seriously frowned upon, so please don't do it. Other recommended usage parameters: For any series of more than 100 requests, please do this on the weekend or outside peak times in the USA.  

**Additionally:** this script does not take locus orientation into account, so manual verification that all spacers are in the **same orientation** is essential. 

#### Usage

`blast_parser.py -h`

```
Usage: blast_parser.py [OPTIONS]

Options:
    REQUIRED:
    -f  input file | in XML format only

    OPTIONAL:
    -c  clean | BOOLEAN | convert XML BLAST results to .csv, removing 'Empty hits'
    -t  trim | BOOLEAN | removes hits NOT containing the text 'plasmid, phage, bacteriophage, or prophage'
    -p  predict | BOOLEAN | return aligned flanking regions from BLAST hit Accession numbers
    -l  local | BOOLEAN | use when the previous BLAST query was made on a local database (-remote option was not used)
```

### Examples
First, confirm all spacers are in the same orientation, then run a BLAST query using [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download), setting the output format to XML (-outfmt 5).

```
blastn -query spacers_concat.fasta -out matches.xml -evalue 1e-3 -remote -db nt -outfmt 5
```
*outputs: **matches.xml** file*<br/><br/>

**Ex 1.** Convert the XML results to CSV using the **-c** (clean) option.
```
blast_parser.py -c -f matches.xml
```

**Ex 2.** Convert the XML results to CSV using the **-c** (clean) option and **-l** (local) option. Use **-l** when blast was executed against a **local** database (when -remote flag was not used during blast query). This simply adjusts how the blast results are parsed.
```
blast_parser.py -cl -f matches.xml
```

**Ex 3.** Convert XML results to CSV and remove BLAST matches not containing the text 'phage, plasmid, bacteriophage, or prophage' (**-t** option).
```
blast_parser.py -ct -f matches.xml
```

**Ex 4.** Predict the PAM sequence from spacer hits (**-p** option).
```
blast_parser.py -ctp -f matches.xml
```
The **-p(x)** option generates a **spacer alignment file** (pam_predict_spacer_list.fa) containing aligned spacers with their corresponding flanks, and a **flank alignment file** (pam_predict_flanks.fa), in which the spacer sequence has been removed and replaced with --, for easier flank alignment. Open these files in an alignment visualization tool, like [IGV](http://software.broadinstitute.org/software/igv/) or [Geneious](https://www.geneious.com/).

Flank alignment file example:</br>
<img src="https://github.com/CRISPRlab/CRISPRutils/blob/master/img/PAM.png" width="400">




### Pro Tip
The **CSV** file created using the **-c** option can be **modified manually**. Open the CSV file in your text editor of choice, then delete any rows with BLAST hits not to your liking (organisms you want to exclude, etc.). Make sure to save the file when complete, but do not rename it. Then run blast_parser.py with the **-p(x)** options. It will pickup the CSV file changes and search Entrez accordingly.
