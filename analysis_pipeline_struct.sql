------------------------------------------------------------
-- 1. Outcome: reshape optimal functioning to long format
------------------------------------------------------------
DROP TABLE IF EXISTS outcome_long;

CREATE TABLE outcome_long AS
SELECT
    participant,
    'B0' AS session,
    ws_cpa AS optimal_functioning
FROM acapcnsvs
WHERE ws_cpa IS NOT NULL
UNION ALL
SELECT
    participant,
    'M3',
    ws_c3m
FROM acapcnsvs
WHERE ws_c3m IS NOT NULL
UNION ALL
SELECT
    participant,
    'M6',
    ws_c6m
FROM acapcnsvs
WHERE ws_c6m IS NOT NULL;

------------------------------------------------------------
-- 2. Demographics + time variables (filter NULLs)
------------------------------------------------------------
DROP TABLE IF EXISTS demo_clean;

CREATE TABLE demo_clean AS
SELECT
    participant,
    session,
    dpi_mri AS time,
    dpi_mri * dpi_mri AS time2,
    ageinjury AS age,
    gender AS sex,
    CASE
        WHEN grp = 1 THEN 'Concussion'
        WHEN grp = 2 THEN 'Orthopedic'
        ELSE 'Unknown'
    END AS grp
FROM acap_demo
WHERE dpi_mri IS NOT NULL
  AND session IS NOT NULL;

------------------------------------------------------------
-- 3. Join demo with outcome
------------------------------------------------------------
DROP TABLE IF EXISTS demo_outcome;

CREATE TABLE demo_outcome AS
SELECT
    d.participant,
    d.session,
    d.time,
    d.time2,
    d.age,
    d.sex,
    d.grp,
    o.optimal_functioning
FROM demo_clean d
LEFT JOIN outcome_long o
    ON d.participant = o.participant
   AND d.session = o.session;

------------------------------------------------------------
-- 4. Final analysis table (metrics + demo + outcome)
------------------------------------------------------------
DROP TABLE IF EXISTS analysis_data_struct;

CREATE TABLE analysis_data_struct AS
SELECT
    g.participant,
    g.session,
    d.time,
    d.time2,
    d.age,
    d.sex,
    d.grp,
    g.metric,
    g.metric_value,
    g.roi,
    g.node,
    g.hemisphere,
    d.optimal_functioning
FROM gt_struct g
JOIN demo_outcome d
    ON g.participant = d.participant
   AND g.session = d.session
WHERE g.thresh = 'weighted'
  AND g.metric_value IS NOT NULL;

------------------------------------------------------------
-- 5. Indexes for performance
------------------------------------------------------------
CREATE INDEX IF NOT EXISTS idx_analysis_pid
    ON analysis_data_struct(participant);

CREATE INDEX IF NOT EXISTS idx_analysis_metric
    ON analysis_data_struct(metric);

CREATE INDEX IF NOT EXISTS idx_analysis_session
    ON analysis_data_struct(session);

	
	
	
