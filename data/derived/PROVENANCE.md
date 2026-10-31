# Provenance of the committed derived tables

## climate_worldclim_10m.csv

6,897 georeferenced records with WorldClim 2.1 BIO1–BIO19 at 10 arc-minute
resolution, and a collection year where one could be recovered (4,736 records,
69%).

This table was carried over from the preceding reference-library reliability
study rather than regenerated here, so that the analysis runs without a fresh
set of downloads. Its occurrence component is a GBIF `occurrence/search` pull
for *Tetranychus urticae* alone, which carries **no citable download DOI**.

Before the manuscript is submitted, replace it:

```
python scripts/fetch_gbif.py "Tetranychus urticae" gbif_tetranychus_urticae --download
python scripts/fetch_gbif.py "Tetranychus" gbif_tetranychus_genus --download
python scripts/fetch_gbif.py "Oligonychus" gbif_oligonychus_genus --download
python scripts/clean_coords.py
python scripts/climate_extract.py 10m
```

That version will carry a GBIF download DOI, will have passed
`clean_coords.py`, and will include the second genus. Record counts and every
number in `results/` will change.

| Column group | Source |
|---|---|
| `source`, `record_id`, `taxon`, `country`, `host`, `year`, `lat`, `lon` | NCBI GenBank, BOLD, GBIF |
| `bio1`–`bio19` | WorldClim 2.1, 10 arc-minute (Fick & Hijmans 2017) |

Sequence-record coordinates derived from country or province centroids were
excluded at extraction. GBIF coordinates in this version have **not** yet been
through `clean_coords.py`.
