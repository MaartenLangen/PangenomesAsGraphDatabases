import argparse
import pandas as pd
import os
import sys 
import csv
from datetime import datetime

def getFeatureDf(pirateDir: str) -> pd.DataFrame:
    """Get all relevant info from features, except for the sequences.

    ## Input
    - pirateDir: string refering to directory with PIRATE output.

    ## Output
    - featureDf: a pd.DataFrame with Name, Start, End, Length, Strand, 
    Product, Strain, and FeatureType for each feature
    """

    # Get all relevant info from co-ords folder
    featureDf = pd.DataFrame()
    coordsDir = os.path.join(pirateDir, "co-ords")
    filenames = os.listdir(coordsDir)
    for filename in filenames:
        file = os.path.join(coordsDir, filename)
        newCoordsDf = pd.read_csv(file, sep="\t", header=0, usecols=["Name", "Start", "End", "Length", "Type", "Strand", "Product"])
        newCoordsDf["Strain"] = pd.Series([filename.split(".")[0] for x in range(len(newCoordsDf.index))])
        featureDf = pd.concat([featureDf, newCoordsDf])
    featureDf.rename(columns= {'Type': 'FeatureType'}, inplace=True)
    return featureDf

def getReprSeqDf(pirateDir: str) -> pd.DataFrame:
    """ Get a dataframe which for each cluster saves the representative locus.
    
    ## Input: 
    - pirateDir: string refering to directory with PIRATE output.
    
    ## Output
    - reprSeqDf: pd.DataFrame with cluster_name and reference_locus
    """

    # Get representative loci, as defined by PIRATE from representative_sequences.ffn
    file = open(os.path.join(pirateDir, "representative_sequences.ffn"))
    representativeSequences = {"cluster_name": [], "reference_locus": []}
    for line in file:
        if line.count(">") > 0: # indicates heading in fasta file
            cluster_name = line.split(";")[0].replace(">","")
            representativeSequences["cluster_name"].append(cluster_name)
            reference_locus = line.split(';')[2].replace('locus_tag=','')
            representativeSequences["reference_locus"].append(reference_locus)
    file.close()
    # convert the dictionary to a pd.DataFrame
    reprSeqDf = pd.DataFrame(representativeSequences)

    return reprSeqDf

def getClusterDf(pirateDir: str, reprSeqDf: pd.DataFrame) -> pd.DataFrame:
    """ Get all relevant info from clusters, except for the reference sequences.

    ## Input
    - pirateDir: string refering to directory with PIRATE output.
    - reprSeqDf: pd.DataFrame with cluster_name and reference_locus.

    ## Output
    - clusterDf: pd.DataFrame with allele_name, gene_family (cluster), concsensus_product,
    threshold, number_genomes, min_length, max_length, average_length, feature, and reference_locus
    """

    # location of file with cluster info
    file = os.path.join(pirateDir, "PIRATE.gene_families.tsv")
    # Determine the total nb of columns in this file:
    num_cols = 0
    with open(file) as f:
        reader = csv.reader(f, delimiter='\t', skipinitialspace=True)
        first_row = next(reader)
        num_cols = len(first_row)
    # relevant columns to read
    cols = [0,1,3,4,6]
    cols.extend(range(17,num_cols))
    # read file as pd.DataFrame
    clusterDf = pd.read_csv(file, sep="\t", header=0, usecols=cols)
    # Convert column with genes per strain to one column with all genes
    for col in range(9,clusterDf.shape[1]):
        clusterDf.iloc[:,8] = clusterDf.iloc[:,8].fillna("0") + ";" + clusterDf.iloc[:,col].fillna("0")
    # rename some columns to indicate their new meaning and make them readable in Neo4j
    clusterDf.drop(clusterDf.columns.values.tolist()[9:num_cols], axis=1, inplace=True)
    clusterDf.rename(columns = {"min_length(bp)": "min_length", "max_length(bp)": "max_length", "average_length(bp)": "average_length", clusterDf.columns[8]: "feature"}, inplace = True)
    # Format indicating that strain has multiple features matching to cluster should be adjusted to make cypher-readable
    clusterDf['feature'] = clusterDf['feature'].apply(lambda x: x.replace("(", "").replace(")", "").replace(":",";"))
   
    # save representative locus to the cluster dataframe
    clusterDf = clusterDf.merge(reprSeqDf, how = 'inner', left_on = 'gene_family', right_on = 'cluster_name', validate = 'many_to_one')
    # drop this column, as it holds same info as 'gene_family'
    clusterDf.drop('cluster_name', axis = 1, inplace=True)

    return clusterDf

def getAlignedSeqDf(pirateDir: str, reprSeqDf: pd.DataFrame) -> pd.DataFrame:
    """Get sequences for all the genes aligned to reference 

    ## Input
    - pirateDir: string refering to directory with PIRATE output.
    - reprSeqDf: pd.DataFrame with cluster_name and reference_locus.

    ## Output
    - alignedSeqDf: pd.DataFrame with Feature, Seq (aligned sequence), Cluster, 
    reference_locus, and Seq_reference (sequence of the reference locus)
    """

    # Get sequences for all the features aligned to reference 
    alignedSeq = {"Feature": [], "Seq": [], "Cluster": []}
    seqDir = os.path.join(pirateDir, "feature_sequences")
    filenames = os.listdir(seqDir)

    for filename in filenames:
        if filename.split(".")[1] != "nucleotide": continue
        file = os.path.join(seqDir, filename)
        fasta = open(file)
        for line in fasta:
            if line.count(">") > 0: 
                alignedSeq["Feature"].append(line.replace(">","").replace("\n",""))
                alignedSeq["Cluster"].append(filename.split('.')[0])
            else: 
                alignedSeq['Seq'].append(line.replace("\n", ""))
        fasta.close()

    alignedSeqDf = pd.DataFrame.from_dict(alignedSeq)

    alignedSeqDf = alignedSeqDf.merge(reprSeqDf[["cluster_name", 'reference_locus']], how='inner', left_on='Cluster', right_on='cluster_name', validate='many_to_one')

    # add column with reference sequence of the relevant cluster
    alignedSeqDf = alignedSeqDf.merge(alignedSeqDf[["Feature", "Seq"]], how='inner', left_on='reference_locus', right_on='Feature', suffixes=["", "_reference"])
    alignedSeqDf = alignedSeqDf.drop(["cluster_name", "Feature_reference"], axis=1)

    return alignedSeqDf

def getStrainSeqDf(pirateDir: str) -> pd.DataFrame:
    ''' Get a dataframe with full sequences per strain.

    ## Input
    - pirateDir: string refering to directory with PIRATE output.

    ## Output
    - strainSeqDf: a pd.DataFrame with Strain and Sequence.
    '''

    strainSequences = {'Strain': [], 'Sequence': []}
    seqDir = os.path.join(pirateDir, "modified_gffs")
    filenames = os.listdir(seqDir)

    for filename in filenames:
        file = os.path.join(seqDir, filename)
        sequence = ''
        fasta = open(file)
        fastaStarted = False
        for line in fasta:
            # many lines with extra info can be skipped here
            if not fastaStarted: 
                # marks the beginning of the sequence
                if line.count('##FASTA') > 0: 
                    fastaStarted = True
            else:
                # contigs are seperated by >...
                if not line.count(">") > 0:
                    sequence = sequence + line
        fasta.close()
        strainSequences['Strain'] += [filename.split('.')[0]]
        strainSequences['Sequence'] += [sequence.replace('\n','')]

    strainSeqDf = pd.DataFrame.from_dict(strainSequences)
    # this index makes accessing sequence for specific strain easier
    strainSeqDf.index = strainSeqDf['Strain']

    return strainSeqDf

def determineVariation(seq: str, refSeq: str) -> str:
    """ Compares seq to refSeq and returns the differences as a string.

    ## Input
    - seq: str with sequence to be compared
    - refSeq: str with sequence to compare to

    ## Output
    - DELTA: 'comressed' difference. Format: alternating number and nucleotide,
    with number indicating distance from previous difference.
    """

    previousDiff = 0
    DELTA = ''
    if seq is refSeq: 
        return DELTA
    else: 
        for i in range(0,len(seq)):
            if seq[i] != refSeq[i]:
                dist = i - previousDiff
                DELTA = DELTA + str(dist) + seq[i]
                previousDiff = i
        return DELTA

def addVariationToFeatureDf( alignedSeqDf: pd.DataFrame, featureDf: pd.DataFrame) -> pd.DataFrame:
    """ Adds column to featureDf with difference to reference Seq (called variation)
    
    ## Input
    - alignedSeqDf: pd.DataFrame with feature, Seq (aligned sequence), Cluster, 
    reference_locus, and Seq_reference (sequence of the reference locus)
    - featureDf: a pd.DataFrame with Name, Start, End, Length, Strand, 
    Product, Strain, and FeatureType for each feature
    
    ## Output
    - featureType: an updated pd.DataFrame with Name, Start, End, Length, Strand, 
    Product, Strain, FeatureType, and Variation for each feature
    """

    # create empty series that will hold variation strings
    variationSeries = pd.Series("", alignedSeqDf.index, name='Variation')
    # iterate over clusters
    for cluster in alignedSeqDf["Cluster"].unique():
        # this dictionary will map sequences to their variation as compared 
        # to reference sequence of the cluster
        variations = {}
        # rows related to current cluster
        rows = alignedSeqDf.loc[alignedSeqDf["Cluster"] == cluster].index
        # the reference sequence of this cluster
        refSeq = alignedSeqDf.loc[rows[0], "Seq_reference"]
        # list of unique sequences in this cluster
        sequences = alignedSeqDf.loc[rows,"Seq"].unique()
        # iterate over unique sequences in cluster
        for seq in sequences:
            # add unique sequence as key to disctionary and map to variation
            variations[seq] = determineVariation(seq, refSeq)
        # save the variations for each feature to the series
        variationSeries.iloc[rows] = alignedSeqDf.loc[rows,"Seq"].map(variations)
    # merge this variation in alignedSeqDf (has to be done to make sure that
    # each feature had the right variation, as these objects have the same indexes)
    alignedSeqDf = pd.merge(alignedSeqDf['Feature'], variationSeries, how = 'inner', left_index = True, right_index = True)
    # merge alignedSeqDf in featureDf so that it has information on variation. Left
    # merge is performed, as there are some sequences without clusters and thus
    # without alignedSeqDfs.
    featureDf = featureDf.merge(alignedSeqDf, how = 'left', left_on = 'Name', right_on = 'Feature', validate='one_to_one', indicator='featureDf.merge(alignedSeqDf)')
    # remove this column, as it is redundant and only used for merging
    featureDf = featureDf.drop('Feature', axis=1)

    return featureDf

def addFullSequencesToFeatureDf(featureDf: pd.DataFrame, strainSeqDf: pd.DataFrame) -> pd.DataFrame:
    '''Some sequences are not part of any cluster. Add their full sequence in 
    a seperate column.

    ## Input
    - featureDf: a pd.DataFrame with Name, Start, End, Length, Strand, 
    Product, Strain, FeatureType, and Variation for each feature.
    - strainSeqDf: a pd.DataFrame with Strain and Sequence.

    ## Output
    - featureDf: an updated pd.DataFrame with Name, Start, End, Length, Strand, 
    Product, Strain, FeatureType, Variation and FullSequence for each feature (either Variation
    or FullSequence hold info on sequence).
    '''

    # Series where feature sequences will be stored until merged to featureDf
    fullSequenceSerie = pd.Series('', featureDf.index, name='FullSequences')
    # get indexes of rows without variation info
    rows = featureDf.loc[featureDf['featureDf.merge(alignedSeqDf)'] == 'left_only'].index
    # iterate over these rows and add sequences
    for rowIndex in rows:
        strain = featureDf.loc[rowIndex, 'Strain']
        start = featureDf.loc[rowIndex, 'Start']
        end = featureDf.loc[rowIndex, 'End']
        fullSequence = strainSeqDf.loc[strain, 'Sequence']
        fullSequenceSerie.iloc[rowIndex] = fullSequence[start-1:end]
    
    # merge the sequences to the featureDf and delete column indicating whether sequence was missing
    featureDf = pd.merge(featureDf,fullSequenceSerie, how='left', left_index=True, right_index=True)
    featureDf = featureDf.drop('featureDf.merge(alignedSeqDf)', axis=1)

    return featureDf

def addReferenceSeqToClusterDf(clusterDf: pd.DataFrame, alignedSeqDf: pd.DataFrame) -> pd.DataFrame:
    """ Add reference sequences to clusterDf.

    ## Input
    - clusterDf: pd.DataFrame with allele_name, gene_family (cluster),
    concsensus_product, threshold, number_genomes, min_length, max_length,
    average_length, feature, and reference_locus
    - alignedSeqDf: pd.DataFrame with Feature, Seq (aligned sequence), Cluster, 
    reference_locus, and Seq_reference (sequence of the reference locus)

    ## Output
    - clusterDf: updated pd.DataFrame with allele_name, gene_family (cluster),
    concsensus_product, threshold, number_genomes, min_length, max_length,
    average_length, feature, reference_locus, and Seq_reference
    """

    clusterDf = clusterDf.merge(alignedSeqDf[['Feature', 'Seq']], how='left', left_on = 'reference_locus', right_on = 'Feature', suffixes = ['','_reference'])
    # remove this column, as it is redundant and only used for merging
    clusterDf.drop(['Feature'], axis = 1, inplace = True)

    return clusterDf

def getNeighbourEdgesDf(featureDf: pd.DataFrame) -> pd.DataFrame:
    """Use start positions of genes to order them and determine neighbours
    
    ## Input
    - featureDf: an updated pd.DataFrame with Name, Start, End, Length, Strand, 
    Product, Strain, FeatureType, Variation and FullSequence for each feature.

    ## Output
    - neighbourEdgesDf: a pd.DataFrame with source and receiving feature
    """

    # sort relevant columns of featureDf
    sortedFeatureDf = featureDf[['Name','Strain', 'Start']].sort_values(by=['Strain','Start'])
    # initiate some variables needed across for-loop iterations
    previous_id = ""
    previous_start = 10001 # large enough so that first gene is has smaller start
    neighborEdgesFeature = {'sourceFeature': [], 'receivingFeature': []} # dict will contain origin gene id and destination gene id of neighbor edges
    for row in sortedFeatureDf.iterrows():
        [name, strain, start] = list(row[1])
        if previous_start < start: # indicates that we're looking at same strain -> so neighbors
            neighborEdgesFeature['sourceFeature'].append(previous_id)
            neighborEdgesFeature['receivingFeature'].append(name)
        # set parameters for next iteration
        previous_id = name
        previous_start = start
    # convert dict to dataframe
    neighbourEdgesDf = pd.DataFrame(neighborEdgesFeature)

    return neighbourEdgesDf

def exportDfToCSV(df: pd.DataFrame, file: str):
    df.to_csv(file, sep=",", index=False)
    

def __main__():
    # parse arguments
    parser = argparse.ArgumentParser(description='Convert pirate output to Neo4j readable CSV files.')
    parser.add_argument('pirateDir', type=str, help='Directory containing pirate output')
    parser.add_argument('outputDir', type=str, help='Directory where csv-files should be saved')
    args = parser.parse_args()

    # print start
    print("Start Time = ", datetime.now().strftime("%H:%M:%S"))

    # get initial information from files
    featureDf = getFeatureDf(args.pirateDir)
    print("featureDf made, intitial.  ", datetime.now().strftime("%H:%M:%S"))
    reprSeqDf = getReprSeqDf(args.pirateDir)
    print("reprSeqDf, initial")
    clusterDf = getClusterDf(args.pirateDir, reprSeqDf)
    print("clusterDf, initial")
    alignedSeqDf = getAlignedSeqDf(args.pirateDir, reprSeqDf)
    print("alignedSeqDf,initial")
    strainSeqDf = getStrainSeqDf(args.pirateDir)
    print("strainSeqDf, initial   ", datetime.now().strftime("%H:%M:%S"))

    # Add sequence information to featureDf
    featureDf = addVariationToFeatureDf(alignedSeqDf, featureDf)
    featureDf = addFullSequencesToFeatureDf(featureDf, strainSeqDf)
    print("featureDf have sequence info   ", datetime.now().strftime("%H:%M:%S"))

    # Add sequence information to clusterDf
    clusterDf = addReferenceSeqToClusterDf(clusterDf, alignedSeqDf)
    print("clusterDf have sequence info ", datetime.now().strftime("%H:%M:%S"))

    # Get info on neighbour edges:
    neighbourEdgesDf = getNeighbourEdgesDf(featureDf)
    print("neighbourEdgesDf, ready to export    ", datetime.now().strftime("%H:%M:%S"))

    # export these dataframes
    exportDfToCSV(featureDf, os.path.join(args.outputDir, "featureNodes.csv"))
    exportDfToCSV(clusterDf, os.path.join(args.outputDir, "clusterNodes.csv"))
    exportDfToCSV(neighbourEdgesDf, os.path.join(args.outputDir, "neighbourEdges.csv"))

    # print start
    print("End Time = ", datetime.now().strftime("%H:%M:%S"))


__main__()