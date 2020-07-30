# SNP Primer Design Pipeline 2
These scripts make a simple pipeline to design **KASP** (Kompetitive Allele Specific PCR) and **CAPS/dCAPS** primers for SNP genotyping. Compared to **[SNP Primer Design Pipeline](https://github.com/pinbo/SNP_Primer_Pipeline)** (1), this version should be able to design primers for KASP/CAPS primers for **ANY species**, not just polyploid wheat. I have not tested in other species, so please open a new issue!

Please also send me an email if it works for your species. I will add it here, so others can try without problems. Thanks!


# Usage

I divided the pipeline into 7 steps:
- Script "parse_polymarker_input.py": parse the polymarker input and prepare a fasta file for blast
- Blast using system command "blastn"
- Script "getflanking.py": Parse the blast output file and output the homelog contigs and flanking ranges
- Split the range file for each marker with system command "awk"
- Get flanking sequences for each file with command "blastdbcmd"
- Get KASP primers using script "getkasp3.py"
- Get CAPS primers using script "getCAPS.py"


You can run this step by step or run the whole pipeline with script "run_getkasp.py". I suggest run the 7 steps in the script "run_getkasp.py" step by step to get familiar how it works first. **You can try it first with the examples in the examples folder**. 

Example: `run_getkasp.py for_polymarker.csv 200 1 1 63 25 0 /home/junli/blastdb/iwgsc_refseqv1.0.fasta`

**Inputs are**: 1) SNP files with the polymarker format (each line is a SNP with 3 comma separated fields: SNP name, chromosome/contig where it is located, flanking sequences), 2) enzyme maximum price (per 1000 U), 3) whether to design CAPS (1 for yes and 0 for NO), 4) whether to design KASP (1 for yes and 0 for NO),  5) maximum Tm (63 C for example), 6) maximum primer size (25 bp for example), 7) whether to pick primer anyway even if it violates specific constrains (1 for yes and 0 for NO), 8) reference file path.

The "bin" folder has all the scripts for each step and software **primer3** and **muscle** in case your system does not have them.

# SNP file input example

This software uses the same format as polymarker (http://www.polymarker.info/about), i.e. a csv file with each line looks like this:

IWA7892,chr1B,AGGATTCACGGGAAAAGATTTCGTCGCA[C/T]GTGCTAGGGTCTTTGCAGAATGTTA

NO space is allowed. It includes Gene ID, chromosome/contig where it is located, flanking sequence (50 to 100 bp) with the alleles in the middle.

# How it works
1. Find all the different sites that can differ all other sequences from the user provided alignment file;
2. Use these sites and the SNP site as SEQUENCE_FORCE_RIGHT_END in primer3 to design all possible left and right primers in the target sequence.

# Pseudo code
1. Read the polymarker input and get:
	- snp position
	- sequences and make them a fasta file for blast later
2. Blast using the fasta file and output blast results
3. Process blast output file to extract flanking sequences (250 bp each side)
4. Multiple sequence alignment of the homeologs
5. Use the msa file to design primers using primer3

# Main Changes
- 07/23/2020: Now supports guessing the best chromosome/contig location based on the blast result: if the chromosome name in the SNP input file did not have any matches in the blast result, it will use the best hit as target chromosome.
- 07/22/2020: Update all code to run with python3 (use `2to3` to convert python2 script to python3 script)
- 07/22/2020: Update to **SNP Primer Design Pipeline 2** for ANY species (before it was just for wheat. Now should be able to work for any species with any ploidy)
- 10/23/2019: Change all script to use python2 and make it easy for users to implement on there Galaxy server.
- 09/02/2019: add "PRIMER\_PICK\_ANYWAY" option for the situation when no primers were obtained.
- 09/01/2019: add a primer3 global setting file for easy change some parameters.
- 09/01/2019: add maximum primer length parameter for low GC content region
- 08/02/2019: add a maximum hits filter: if more than 6 hits, do not design for this SNP, because some SNPs have too many hits.
- 05/26/2019: updated SNP position to polymarker input to fit BLAST+ 2.9.0+. So you need to update your BLAST+ from 2.6.0 to 2.9.0, because the blastdbcmd output format changed in 2.9.0.
- 03/30/2018 Add a new script getCAPS-with-user-input to use user provided multipe sequences to design CAPS/dCAPS primers and the variations can also be an indel now.

# Dependencies

SNP_Primer_Pipeline needs following 3 software to find differences among homeologs and design primer.
1. **Python 3**: Please install python 3.x (default download will be python 3, because python 2 is not supported anymore)
2. **Muscle**: Multiple sequence alignment program (http://www.drive5.com/muscle/)
3. **Primer3**: program for designing PCR primers (http://primer3.sourceforge.net/)
4. **BLAST+ 2.9.0** (or later) package from NCBI (https://blast.ncbi.nlm.nih.gov/Blast.cgi)

"**muscle**" and "**primer3_core**" are included in the package, so "**BLAST+**" and "**Python 3**" software you need to install in your system.

# How to implement it to your own Galaxy server

**Before doing this, please make your blast database ready.** For wheat, please download the RefSeqv1.0 for example, and extract the fasta file, then run the command below to make it blastable:

`makeblastdb -in 161010_Chinese_Spring_v1.0_pseudomolecules.fasta -dbtype nucl -parse_seqids`

I suggest put both the "**SNP Position to polymarker input**" and the "**SNP Design Pipeline 2**" tools in the Galaxy tool menu.

1. Go to the Galaxy root folder and go to the "tools" folder: `cd tools`

1. Clone the SNP Primer Design Pipeline to the "tools" folder: `git clone https://github.com/pinbo/SNP_Primer_Pipeline2.git`

1. Go to the "SNP\_Primer\_Pipeline" folder: `cd SNP_Primer_Pipeline2`

1. Make a copy of the file "*SNP2polymarker.xml.example*" and rename it "*SNP2polymarker.xml*". Or `cp SNP2polymarker.xml.example SNP2polymarker.xml`. Do the same thing for "*getkasp.xml.example*": `cp getkasp.xml.example getkasp.xml`

1. Edit the configuration files "SNP2polymarker.xml" and "getkasp.xml". At least change the reference file location (Red rectangle in the screenshot below): "value=" is the location.
![change reference location](./files/change-references.png)

1. Go back to the Galaxy root folder and go to the "config" folder and add the tool xml location in the file "*tool\_conf.xml*". If the file is not there, just make a copy of "*tool\_conf.xml.sample*" and rename it "*tool\_conf.xml*". Then add the tool xml location there (**The image below is from version 1, so change "SNP_Primer_Pipeline" to "SNP_Primer_Pipeline2"**):
```{xml}
<section name="SNP_Primer_Design" id="snp_primer">
    <tool file="SNP_Primer_Pipeline2/getkasp.xml" />
    <tool file="SNP_Primer_Pipeline2/SNP2polymarker.xml" />
</section>
```
![add xml file to tool_config](./files/add-xml.png)

Now the tools should be there. 
![tools on galaxy](./files/tools-on-galaxy.png)

More details about adding your own tools to Galaxy can be found here: https://galaxyproject.org/admin/tools/add-tool-tutorial/

# Acknowledgements
I borrowed ideas from the polymarker scripts (https://github.com/TGAC/bioruby-polyploid-tools), a great tool for Genome Specific KASPar design in polyploid species. Thanks to the author of Polymarker.

I also borrowed some codes from biopython (https://github.com/biopython/biopython/blob/master/Bio/Emboss/Primer3.py). Thanks to them too.

Thanks to the open source software **Primer3** (http://primer3.sourceforge.net/), **Muscle** (http://www.drive5.com/muscle/),  blast+ package from NCBI (https://blast.ncbi.nlm.nih.gov/Blast.cgi), and **SNP2CAPS** (http://pgrc.ipk-gatersleben.de/snp2caps/).
