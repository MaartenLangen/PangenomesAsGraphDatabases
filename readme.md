This extended data set contains 770 genomes in gbk format downloaded from NCBI. This website was used to retrieve the data: https://www.ncbi.nlm.nih.gov/datasets/genomes/. Here, the term "Pseudomonas aeruginosa" was used in the search bar. In the filters on the subsequent results page, annotated was switched on and assembly level was set to only include complete genomes. This resulted in the 770 strains listed in `strains.txt`.

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

Pirate was performed with these files (on the VSC-HPC). The following command was used (this was done using the environment shared in `pirateEnv.yml`)

```bash
PIRATE --input ./gff3/ --output ./pirateOutput --nucl --align --features "CDS,tRNA,pseudogene" --threads 18
```

There were some problems with sanity checks. This problem was solved by removing two lines in one of the scripts used by PIRATE (`pangenome_construction.pl`): removed lines 653 and 1005.

The pirate output can be converted to csv-files that are fit to import in a Neo4j database by executing the following command in the folder of this document: (make sure that the second folder in the command already exists, otherwise, there will be an error that the file can't be written) (python code should be executed in conda environment from `pythonEnv.yml`)

```
python pirateToDatabase.py ./02_pirateOutput ./03_graphDBInput
```

Now, an empty database in Neo4j can be created. Remember the password used, as it will be needed later on to allow python scripts to interact with the database (I used '770strains', so look for this string in python scripts with connection to Neo4j and replace it with your own password). Once created, click the three dots in the line of the new database, select ``open folder > Import`` and copy the files from ``03_graphDBInput`` to that location. Now the database can be started. Open a neo4j browser window and execute the cypher code from `02_CypherScripts/01_creatingDatabase.cypher`.

Further annotation can be added to this database for the detection of genomic islands. Clusters with sequence similarity to virus database are annotated using `extendingDB_performingBlast.py` (the graph database should be running and the uri+password should be changed in the script, and the script should be executed in a linux environment with BLAST added to path). Nucleotide statistics for the features and strain averages are calculated in `extendingDB_FeaturesAndStrainMetrics.py`, which generates a file that has to be imported in the neo4j in a similar way as the original data (the graph database should be running and the uri+password should be changed in the script). Cypher code to integrate the data in the graph database is provided in `02_CypherScripts/02_addNucleotideCompositionMetrics.cypher`. 

The analyses are present in `03_analyses`. Again, the graph database should be running and the correct uri+password should be used to allow interaction with the database.