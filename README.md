# reference-coverage

A workflow for measuring how well a public DNA barcode reference library covers
the geographic and climatic range in which a taxon is actually recorded, and for
ranking the places where a new reference specimen would add the most information.

The workflow differs from existing reference-library gap analyses in one respect:
it does not assume that the species label on a public sequence record is correct.
Each record carries a reliability weight derived from whether its label can be
confirmed against morphologically identified anchor specimens. Coverage is then
measured as reliability-weighted density rather than record count, so a climate
zone represented only by doubtful records is reported as uncovered.

Spider mites (Acari: Tetranychidae) are the demonstration taxon. Nothing in the
pipeline is specific to them: every script takes a taxon name as an argument.

## Inputs

| Source | What it provides | How it arrives |
|---|---|---|
| GBIF | Occurrence records with coordinates and dates — the recorded distribution | `fetch_gbif.py` |
| NCBI GenBank, BOLD | COI sequence records, some with coordinates — the reference library | ported from the preceding study |
| WorldClim 2.1 | BIO1–BIO19 at the coordinates of every georeferenced record | `climate_extract.py` |
| Köppen-Geiger (Beck et al.) | Categorical climate zone, for interpretation | `climate_extract.py --koppen` |
| Zenodo 10.5281/zenodo.23808136 | Per-record identification verdicts from the preceding study | `fetch_zenodo.py` |

Rasters and bulk downloads land in `data/raw/` and `data/external/` and are not
versioned. Derived tables small enough to read are committed in `data/derived/`.

## Pipeline

```
fetch_gbif.py        occurrence records for a taxon
clean_coords.py      drop centroid, zero, marine and institution coordinates
climate_extract.py   attach BIO1-BIO19 and Köppen zone to every georeferenced record
fetch_zenodo.py      per-record identification verdicts from the preceding study
deficit.py           PCA climate space, reliability-weighted deficit surface,
                     geographic deficit map, sampling-priority list
holdout.py           build the priority map from pre-2015 records only, then test
                     it against records deposited afterwards
```

## Reliability weights

`deficit.py` maps each identification verdict to a weight in [0, 1]:

| Verdict | Weight | Meaning |
|---|---|---|
| `AGREES` | 1.00 | label confirmed against a morphologically identified anchor in the same cluster |
| `AGREES_SHARED_CLUSTER` | 0.50 | label consistent, but the cluster holds more than one anchored species |
| `NO_ANCHOR` | 0.25 | no anchored specimen in the cluster; label unverifiable |
| `UNNAMED_SUGGESTED` | 0.25 | record unnamed, one species suggested by clustering |
| `UNNAMED_SEVERAL_POSSIBLE` | 0.10 | record unnamed, several species possible |
| `DISAGREES` | 0.00 | label contradicted by an anchored specimen |

These values are a judgement, not a measurement. `deficit.py --sensitivity`
re-runs the whole analysis under alternative schemes and reports how far the
conclusion moves.

## Reproducing

```
pip install -r requirements.txt
python scripts/fetch_gbif.py "Tetranychus urticae" gbif_tetranychus_urticae
python scripts/clean_coords.py
python scripts/climate_extract.py 10m
python scripts/fetch_zenodo.py
python scripts/deficit.py
python scripts/holdout.py 2015
```

## Licence

Code: MIT. Derived data tables in `data/derived/` and `results/`: CC BY 4.0.
GBIF occurrence data carry the licences of their constituent datasets; the
download DOI recorded in `data/raw/gbif_download.txt` is the citable source.
