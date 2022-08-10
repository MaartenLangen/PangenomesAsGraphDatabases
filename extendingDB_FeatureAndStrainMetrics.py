# from ast import Assert
from __future__ import annotations

import sys
from datetime import datetime 
import neo4j
import numpy as np
import pandas as pd
import re
from Bio.Seq import Seq
from Bio.Seq import MutableSeq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import GC
from Bio.SeqUtils import CodonUsage as CU

sys.path.append('/mnt/c/Users/maart/OneDrive/school/2021-2022/Thesis/project/05_extendedData/04_neo4jAPI')
import neo4jConnection

def variationToFullSequence(variation: str, referenceSequence: str) -> str:
    featureSequence = MutableSeq(referenceSequence)
    if variation is not None:
        diffs = re.findall(r'(\d+)(\D)', variation)
        index = 0
        for diff in diffs:
            index += int(diff[0])
            featureSequence[index] = diff[1]
    return str(featureSequence).replace('-','')

def getFeatureVariation_work(tx: neo4j.Session) -> neo4j.Result:
    result = tx.run("""
        MATCH (f:FEATURE {feature_type: 'CDS'})-[:ORTHOLOG]->(c:CLUSTER)
        RETURN f.feature_id, f.variation, c.reference_sequence
    """)

    return result.values('f.feature_id', 'f.variation', 'c.reference_sequence')

def getFeatureSequenceMetrics(driver: neo4j.Driver) -> list[SeqRecord]:
    # create dictionary that will hold relevant info
    records_dict = {'featureID':[],'GC':[], 'CAI': []}
    # retrieve info from database
    with driver.session() as session:
        result = session.read_transaction(getFeatureVariation_work)
        obj = CU.CodonAdaptationIndex() # created to later determine CAI
        for record in result:
            [featureID, variation, referenceSequence] = record

            featureSequence = variationToFullSequence(variation, referenceSequence)

            # determine GC
            GC_cont = GC(Seq(featureSequence))

            # determine CAI
            if 'N' in featureSequence:
                filteredSeq = ''
                i = 0
                while i <= (len(featureSequence)-3):
                    part = featureSequence[i:i+3]
                    if not ('N' in part):
                        filteredSeq = filteredSeq + part
                    i+=3
            elif len(featureSequence) % 3 != 0: 
                filteredSeq = featureSequence[0:len(featureSequence)-len(featureSequence)%3]
            try:
                CAI = obj.cai_for_gene(featureSequence)
            except: 
                CAI = np.nan
            
            # save values to dictionary 
            records_dict['featureID'].append(featureID)
            records_dict['GC'].append(GC_cont)
            records_dict['CAI'].append(CAI)

    # return a pandas dataframe from the dictionary
    return pd.DataFrame.from_dict(records_dict)

def writeFeatureCompositionMetrics_work(tx: neo4j.Session, featureID: str, GC_cont: float, CAI: float) -> neo4j.Result:
    if not (np.isnan(GC_cont) or np.isnan(CAI)): return tx.run("""
        MERGE (f:FEATURE {feature_id: $featureID})
        ON MATCH SET
            f.GC = $GC_cont,
            f.CAI = $CAI
    """, featureID=featureID, GC_cont=GC_cont, CAI=CAI)
    elif np.isnan(CAI): return tx.run("""
        MERGE (f:FEATURE {feature_id: $featureID})
        ON MATCH SET
            f.GC = $GC_cont
    """, featureID=featureID, GC_cont=GC_cont)
    elif np.isnan(GC_cont): return tx.run("""
        MERGE (f:FEATURE {feature_id: $featureID})
        ON MATCH SET
            f.CAI = $CAI
    """, featureID=featureID, CAI=CAI)

def writeStrainsCompositionMetrics_work(tx: neo4j.Session):
    return tx.run("""
        match (s:STRAIN)
        call {
            with s
            match (f:FEATURE)-[:FEATURE_IN_STRAIN]->(s)
            return 
                avg(f.GC) as avg_GC,
                stDev(f.GC) as stDev_GC,
                avg(f.CAI) as avg_CAI,
                stDev(f.CAI) as stDev_CAI
        }
        set 
            s.avg_GC = avg_GC,
            s.stDev_GC = stDev_GC,
            s.avg_CAI = avg_CAI,
            s.stDev_CAI = stDev_CAI
        """)

def writeFeatureAndStrainCompositionMetrics(driver: neo4j.Driver, featuresCompositionDf: pd.DataFrame) -> None:
    with driver.session() as session:
        for index, row in featuresCompositionDf.iterrows():
            session.write_transaction(writeFeatureCompositionMetrics_work, row['featureID'], row['GC'], row['CAI'])
        session.write_transaction(writeStrainsCompositionMetrics_work)


    
def __main__():
    # print start
    print("Start Time = ", datetime.now().strftime("%H:%M:%S"))

    # neo4j connection
    connection = neo4jConnection.Connection(uri, user, pwd)
    driver = connection.getDriver()
    print('database connection made.    ', datetime.now().strftime("%H:%M:%S"))

    featuresCompositionDf = getFeatureSequenceMetrics(driver)
    featuresCompositionDf.to_csv('featuresCompositionDataframe', index=False)
    print("featureCompositionDF obtained.", datetime.now().strftime("%H:%M:%S"))

    # writeFeatureAndStrainCompositionMetrics(driver, featuresCompositionDf)
    # print("finished writing to db.    ", datetime.now().strftime("%H:%M:%S"))

# TODO: Turn this in parameters of the script.
uri = 'bolt://localhost:7687'
user = 'neo4j'
pwd = '770strains'

__main__()