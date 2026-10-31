# reference-coverage

A workflow for measuring how much of the geographic and climatic range in which
a taxon is recorded is backed by *verifiable* reference sequences, and for
ranking the places where a new reference specimen would add the most.

Existing reference-library gap analyses count records. This one counts only what
a user of the library could actually rely on: each sequence record carries a
weight for whether its species label can be confirmed against a morphologically
identified specimen, and coverage is measured as verifiability-weighted density.
A climate zone represented only by unverifiable records is reported as
uncovered.

The claim is about usable coverage, not about accuracy. Across this anchor table
an unchecked label stands about nine times in ten — but nobody reading the
library can tell which tenth is wrong, so an unverifiable record cannot carry the
weight of a reference. Every run prints that base rate as context, and reports a
floor in which every weight is replaced by the bare chance the label is correct,
so the strongest objection to the method is answered in its own output.

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

## What counts as reference coverage

A georeferenced record contributes to coverage only if it supplies a sequence
from the barcode region. `deficit.py` reads the `marker` column and moves
everything else — COI-3P, other genes, records with no sequence — onto the
occurrence side, where it still marks where the taxon was found. The move is
reported, by marker, in every run. A climate table with no `marker` column
keeps all of its sequence records, on the assumption that the query restricted
them to COI.

## Verifiability weights

`deficit.py` maps each identification verdict to a weight in [0, 1]:

| Verdict | Weight | Meaning |
|---|---|---|
| `AGREES` | 1.00 | label confirmed against a morphologically identified anchor in the same cluster |
| `AGREES_SHARED_CLUSTER` | **1/k** | label consistent, but the cluster holds k anchored species, so it fits any one of them. Derived, not chosen; k is 2 throughout the spider mite data, giving 0.50 |
| `NO_ANCHOR` | 0.25 | no anchored specimen in the cluster; label unverifiable |
| `UNNAMED_SUGGESTED` | 0.25 | record unnamed, one species suggested by clustering |
| `UNNAMED_SEVERAL_POSSIBLE` | 0.10 | record unnamed, several species possible |
| `DISAGREES` | 0.00 | label contradicted by an anchored specimen |

Apart from the derived 1/k, these values state what an unverifiable record is
worth as a reference. They are not estimates of how often such labels are right.
`deficit.py --sensitivity`
re-runs the analysis under alternative schemes, and across a range of coverage
thresholds, and reports how far the conclusion moves. The absolute index
depends on the threshold; the naive-to-weighted inflation does not.

Records with no verdict take `DEFAULT_WEIGHT`. Rather than rest on that one
number, every run reports three bounds: those records trusted in full,
discarded, and the floor in which every weight becomes the observed chance that
a label is merely correct. The base rate behind that floor is estimated only
from clusters that contain an anchored specimen — the better-studied part of the
library — so carrying it over to unanchored records is optimistic rather than
neutral, which is why it is reported as a floor and not used as the working
weight.

## One genus at a time

A sequence of one genus supplies no reference coverage for another, so pooling
genera understates the gap. `deficit.py --genus Tetranychus` and
`holdout.py --sweep --genus Tetranychus` restrict both sides and write to a
subfolder of `results/`. Run each genus separately and report them separately;
the pooled run is kept only as a comparison.

## Hold-out validation

`holdout.py <year>` trains on records up to a cutoff and tests against what came
after. `holdout.py --sweep` repeats it across cutoffs from 2006 to 2020, so the
cutoff is a reported robustness result rather than a choice. Cutoffs whose
training library falls below `MIN_TRAIN_SEQ` sequences are flagged as thin and
should be read as indicative only.

## Reproducing

```
pip install -r requirements.txt
python scripts/fetch_gbif.py "Tetranychus urticae" gbif_tetranychus_urticae
python scripts/clean_coords.py
python scripts/climate_extract.py 10m
python scripts/fetch_zenodo.py
python scripts/deficit.py --sensitivity --genus Tetranychus
python scripts/holdout.py --sweep --genus Tetranychus
```

## Licence

Code: MIT. Derived data tables in `data/derived/` and `results/`: CC BY 4.0.
GBIF occurrence data carry the licences of their constituent datasets; the
download DOI recorded in `data/raw/gbif_download.txt` is the citable source.
