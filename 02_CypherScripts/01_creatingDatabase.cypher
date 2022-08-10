// create Features (113013 ms)
:auto USING PERIODIC COMMIT 5000
LOAD CSV with headers FROM "file:///featureNodes.csv" as row
WITH 
    toString(row.Name) AS id,
    toInteger(row.Start) AS feature_start,
    toInteger(row.End) AS feature_end,
    toInteger(row.Length) AS length,
    toString(row.Strand) AS strand, 
    toString(row.Product) AS name,
    toString(row.Strain) AS strain,
    toString(row.Variation) AS variation,
    toString(row.FullSequences) AS full_sequence,
    toString(row.FeatureType) AS feature_type
CREATE (c:FEATURE{
        feature_id: id,
        feature_start: feature_start,
        feature_end: feature_end,
        length: length,
        strand: strand,
        name: name, 
        strain: strain,
        variation: variation,
        full_sequence: full_sequence,
        feature_type: feature_type
});

// create cluster nodes (6050 ms)
:auto USING PERIODIC COMMIT 5000
LOAD CSV with headers FROM "file:///clusterNodes.csv" as row
WITH
    toString(row.allele_name) AS cluster_id,
    toString(row.consensus_product) AS description,
    toInteger(row.threshold) AS threshold,
    toInteger(row.number_genomes) AS number_genomes,
    toInteger(row.min_length) AS min_length,
    toInteger(row.max_length) AS max_length,
    toFloat(row.average_length) AS average_length,
    [x in split(toString(row.feature), ";") WHERE not (x="0")| x] AS feature_ids,
    toString(row.reference_locus) AS reference_locus,
    toString(row.Seq) AS reference_sequence
CREATE (c:CLUSTER{
    cluster_id: cluster_id,
    description: description,
    threshold: threshold,
    number_genomes: number_genomes,
    min_length: min_length,
    max_length: max_length,
    average_length: average_length,
    feature_ids:  feature_ids,
    lonely_cluster: FALSE,
    reference_locus: reference_locus,
    reference_sequence: reference_sequence
});

// create indexes (39 ms)
create index feature_id_index for (f:FEATURE) on f.feature_id;
create index cluster_id_index for (c:CLUSTER) on c.cluster_id;

// create ortholog relation (99357 ms)
:auto MATCH (c:CLUSTER)
UNWIND c.feature_ids AS ID
WITH c, ID
CALL { 
    WITH c, ID
    MATCH (f:FEATURE{feature_id: ID})
    CREATE (f)-[:ORTHOLOG]->(c)
} IN TRANSACTIONS OF 10000 ROWS;

// create clusters for lonely FEATURES (6642 ms)
:auto MATCH (f:FEATURE)
WHERE NOT ((f)-[:ORTHOLOG]-(:CLUSTER))
WITH f
CALL { 
    WITH f
    CREATE (f)-[:ORTHOLOG]->(c:CLUSTER{
        cluster_id: f.feature_id,
        description: f.name,
        threshold: 50,
        number_genomes: 1,
        min_length: f.length,
        max_length: f.length,
        average_length: f.length,
        CDS_ids:  f.feature_id,
        lonely_cluster: TRUE,
        reference_sequence: f.full_sequence
    })
    SET f.full_sequence = NULL
} IN TRANSACTIONS OF 10000 ROWS;

// create neighbour edges between clusters (898094 ms)
:auto USING PERIODIC COMMIT 5000
LOAD CSV with headers FROM "file:///neighbourEdges.csv" as row
WITH
    toString(row.sourceFeature) AS source,
    toString(row.receivingFeature) AS destination
MATCH 
    (c1:CLUSTER)-[:ORTHOLOG]-(f1:FEATURE{feature_id: source}),
    (c2:CLUSTER)-[:ORTHOLOG]-(f2:FEATURE{feature_id: destination})
MERGE (c1)-[r:NEIGHBOUR]->(c2)
    ON CREATE SET
        r.number_of_members = 1,
        r.members = [f1.strain]
    ON MATCH SET
        r.number_of_members = r.number_of_members + 1,
        r.members = r.members + [f1.strain];
		
// create neighbour edges between features ()
:auto USING PERIODIC COMMIT 10000
LOAD CSV with headers FROM "file:///neighbourEdges.csv" as row
WITH
    toString(row.sourceFeature) AS source,
    toString(row.receivingFeature) AS destination
MATCH 
    (f1:FEATURE{feature_id: source}),
    (f2:FEATURE{feature_id: destination})
MERGE (f1)-[r:NEIGHBOUR]->(f2)
SET
        r.strain = f1.strain;

// set cluster.feature_type (75000 ms)
MATCH (f:FEATURE {feature_type: 'pseudogene'})-[:ORTHOLOG]->(c:CLUSTER)
SET c.feature_type = 'pseudogene';
MATCH (f:FEATURE {feature_type: 'tRNA'})-[:ORTHOLOG]->(c:CLUSTER)
SET c.feature_type = 'tRNA';
MATCH (f:FEATURE {feature_type: 'CDS'})-[:ORTHOLOG]->(c:CLUSTER)
SET c.feature_type = 'CDS';

// clusters with NaN in alignedSeqDb (2226 ms)
match (c:CLUSTER)
where c.reference_sequence is null
match (c)<-[:ORTHOLOG]-(f:FEATURE {feature_id: c.reference_locus})
set c.reference_sequence = f.full_sequence

// unique strain names constraint (215 ms)
create constraint unique_strain_name for (s:STRAIN) require s.name is unique;

// create strain nodes and connect them to features (85502 ms)
:auto match (f:FEATURE)
with f, f.strain as strain
call {
    with f, strain
    merge (s:STRAIN {name: strain})
    merge (f)-[:FEATURE_IN_STRAIN]->(s)
} IN TRANSACTIONS OF 10000 ROWS;

// connect clusters to strain nodes (518406 ms)
:auto match (c:CLUSTER)<-[:ORTHOLOG]-(f:FEATURE)-[:FEATURE_IN_STRAIN]->(s:STRAIN)
with c, s 
call {
    with c, s
    merge (c)-[:CLUSTER_IN_STRAIN]->(s)
} in transactions of 10000 rows;