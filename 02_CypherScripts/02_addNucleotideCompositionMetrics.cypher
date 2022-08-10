// Make sure to first run extendingDB_FeatureAndStrainMetrics.py and place the resulting file in the import folder of the neo4j database

// add nucleotide composition metrics to features (254083 ms)
:auto USING PERIODIC COMMIT 1000
LOAD CSV with headers FROM "file:///featuresCompositionDataframe.csv" as row
WITH    
    toString(row.featureID) AS ID,
    toFloat(row.GC) AS GC,
    toFLoat(row.CAI) AS CAI
MERGE (f:FEATURE {feature_id: ID})
ON MATCH SET f.GC = GC, f.CAI = CAI;

// calculate strain averages and stDev (13555 ms)
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
    s.stDev_CAI = stDev_CAI;