This extended data set contains 770 genomes in gbk format downloaded from NCBI. This website was used to retrieve the data: https://www.ncbi.nlm.nih.gov/datasets/genomes/. Here, the term "Pseudomonas aeruginosa" was used in the search bar. In the filters on the subsequent results page, annotated was switched on and assembly level was set to only include complete genomes. This resulted in the 770 genomes present here.

All gbk files had to be moved in one folder for ease of use further on in the pipeline. This was done by executing the following commands in the downloaded folder:

```
for dir in */ 
do 
	new_file_name=$(echo $dir | awk '{print substr( $0, 1, length($0)-1)}')
	mv $dir*.gbff ../../../01_gbff/$new_file_name.gbk
done
```

All files were zipped:

```
gzip ./*.gbk
```

These zipped gbk files can be converted to gff:

```
perl bp_genbank2gff3.pl --dir 01_gbff --outdir 02_gff3raw
```

All files are renamed from `*.gbk.gz.gff` to `.gff`  by executing the following command in `02_gff3raw`:

```
for file in *; do mv "$file" "${file/.gbk.gz.gff/.gff}"; done
```

All files are then converted to prokka gff by executing the following command in `02_gff3raw`: (Notice that this script is different from the original. It now not only includes CDSs, but also tRNA and pseudogene features. Edits are in line 67 and 73) 

```
for file in *; do python ../convert_refseq_to_prokka_gff.py -i $file -o ../03_gff3/$file; done
```

Pirate was performed with these files (on the VSC-HPC). The following command was used 

```bash
PIRATE --input ./gff3/ --output ./pirateOutput --nucl --align --features "CDS,tRNA,pseudogene" --threads 18
```

There were some problems with sanity checks. This problem was solved by removing two lines in one of the scripts used by python: removed lines 653 and 1005.

The pirate output can be converted to csv-files that are fit to import in a Neo4j database by executing the following command in the folder of this document: (make sure that the second folder in the command already exists, otherwise, there will be an error that the file can't be written)

```
python pirateToDatabase.py ./02_pirateOutput ./03_graphDBInput
```

Now, an empty database in Neo4j can be created. Remember the password used, as it will be needed later on to allow python scripts to interact with the database. Once created, click the three dots in the line of the new database, select ``open folder > Import`` and copy the files from ``03_graphDBInput`` to that location. Now the database can be started. Open a neo4j browser window and execute the following code:


Further annotation can be added to this database for the detection of genomic islands. This is done in the `performingBlast.py`script in folder `05_GenomicIslands`. Here, blast alignment is performed against a database with virus genomes and a database with mobile genetic elements. With this script, sequence composition data is also added to the database. While the database is running, perform the following code (after changing the uri and password to correspond to the running database) in `05_GenomicIslands`: (do this in an environment where all packages are present)

```
export PATH=ncbi-blast-2.12.0+/bin/:$PATH
./performingBlast.py
```