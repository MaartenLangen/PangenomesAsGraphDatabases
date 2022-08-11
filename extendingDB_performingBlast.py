import sys
from typing import Union

from datetime import datetime
import neo4j
import pandas as pd
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

sys.path.append('./04_neo4jAPI')
import neo4jConnection

uri = 'bolt://localhost:7687' # make sure correct uri
user = 'neo4j'
pwd = '770strains' # make sure correct password


def getClusterSequences_work(tx: neo4j.Session) -> neo4j.Result:
    """Fuction performs a query on the database to retrieve clusterIDs and 
    reference sequences of all the clusters.

    ## Input
    - neo4j.Session

    ## Output
    - List with list per cluster. The nested list has info on cluster_id and
    reference_sequence.
    """

    result =  tx.run("""
        MATCH (c:CLUSTER)
        WHERE c.reference_sequence IS NOT NULL
        RETURN c.cluster_id, c.reference_sequence, c.description, c.feature_type
    """)

    return result.values('c.cluster_id', 'c.reference_sequence', 'c.description', 'c.feature_type')

def getClusterSequences(driver: neo4j.Driver) -> list[SeqRecord]:
    """Function starts a session and retrieves cluster_ids and 
    reference_sequences for all cluster nodes.

    ## Input
    - driver: a neo4j.Driver holding connection details to the relevant database.

    ## Output
    - my_records: list of SeqRecord (containing Seq, id, description, annotation)
    """
    my_records = []
    with driver.session() as session:
        result = session.read_transaction(getClusterSequences_work)
        for record in result:
            [raw_clusterId, raw_seq, raw_description, raw_feature_type] = record
            seqRecord = SeqRecord(Seq(raw_seq.replace('-','')), id=raw_clusterId, description=raw_description,
               annotations={'feature_type': raw_feature_type})
            my_records.append(seqRecord)
    return my_records
        


def runBlast(query: str, filename_output: str) -> None:
    """ Run BLAST on a query file using an available database and save output 
    in xml-file.

    ## Input
    - query: str referencing fasta-file containingsequences to be queried with 
    BLAST.
    - filename_output: str that will be the filename for the xml-file that is 
    the result of BLAST.

    ## Post
    - after executing this function succesfully, there will be a file 
    filename_output that holds the results of the BLAST run.
    """
    cline = NcbiblastnCommandline(query=query, db='ref_viruses_rep_genomes',
                                  out=filename_output, outfmt=5, num_threads=12)
    stdout, stderr = cline()
    if stdout != '' or stderr != '': raise AssertionError('Executing BLAST should not lead to results', stdout, stderr)

def readBlast(filename: str) -> pd.DataFrame:
    """ Read the xml-file that was the output of a BLAST run and save the relevant information.

    ## Input
    - filename: str that will be the filename for the xml-file that is 
    the result of BLAST.

    ## Output
    - blast_descriptions_df: pd.DataFrame with info on the description of each BLAST hit.
    """
    result_handle = open(filename)
    blast_records = NCBIXML.parse(result_handle)
    blast_descriptions = {'cluster_id': [], 'title': [], 'score': [], 'e': [], 'num_alignments': []}
    for blast_record in blast_records:
        cluster_id = blast_record.query.split(' ')[0] # Allows blast results to be merged with input queries
        for description in blast_record.descriptions:
            blast_descriptions['cluster_id'].append(cluster_id)
            blast_descriptions['title'].append(description.title)
            blast_descriptions['score'].append(description.score)
            blast_descriptions['e'].append(description.e)
            blast_descriptions['num_alignments'].append(description.num_alignments)
    blast_descriptions_df = pd.DataFrame.from_dict(blast_descriptions)

    return blast_descriptions_df

def writeBlastResultVirus_work(tx: neo4j.Session, cluster_id: str, blastTitles, blastScores) -> neo4j.Result:
    """
    Function can be used to write the results of BLAST against the virus 
    database to the graph database. Relevant clusters will be updated to 
    contain properties hasBlastVirusHit (Boolean), blastTitleVirus (list of 
    titles), blastScoresVirus (list of scores).

    ## Input
    - neo4j.Session: a session object
    - cluster_id: unique identifier (string) of the cluster with BLAST hits
    - blastTitles: list of titles returned by BLAST (string)
    - blastScores: list of scores returned by BLAST (integer/float)
    """
    return tx.run("""
        MATCH (c:CLUSTER {cluster_id: $cluster_id})
        SET
            c.hasBlastVirusHit = TRUE, 
            c.blastTitlesVirus = $blastTitles,
            c.blastScoresVirus = $blastScores
    """, cluster_id = cluster_id, blastTitles = blastTitles, blastScores = blastScores) 

def writeBlastResultToDatabase(blast_description: pd.DataFrame, driver: neo4j.Driver) -> None:
    """
    Wrapper function that opens a neo4j session with the supplied driver and
    writes a series of BLAST results to the graph database.

    ## Input
    - blast_description: pd.DataFrame returned by readBlast
    - driver: a neo4j driver object
    """

    with driver.session() as session:
        for cluster in blast_description['cluster_id'].unique():
            blastTitles = blast_description[blast_description['cluster_id']==cluster]['title'].tolist()
            blastScores = blast_description[blast_description['cluster_id']==cluster]['score'].tolist()
            session.write_transaction(writeBlastResultVirus_work, cluster, blastTitles, blastScores)


    
def __main__():

    # print start
    print("Start Time = ", datetime.now().strftime("%H:%M:%S"))

    # neo4j connection
    connection = neo4jConnection.Connection(uri, user, pwd)
    driver = connection.getDriver()
    print('database connection made. ', datetime.now().strftime("%H:%M:%S"))

    # get cluster sequences from database
    my_records = getClusterSequences(driver)
    SeqIO.write(my_records, 'cluster_sequences.fasta', 'fasta')
    print('obtained records.    ', datetime.now().strftime("%H:%M:%S"))
    
    # run blast and save results in database
    try:
        runBlast('cluster_sequences.fasta', 'blastResults.xml')
    except AssertionError as err: print(err.args)
    blast_description = readBlast('blastResults.xml')
    print('blast results virus obtained.    ', datetime.now().strftime("%H:%M:%S"))
    writeBlastResultToDatabase(blast_description, driver)
    print('written virus hits to db.    ', datetime.now().strftime("%H:%M:%S"))

__main__()