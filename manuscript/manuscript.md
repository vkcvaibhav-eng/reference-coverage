# Verifiable coverage of public COI reference libraries in climate space: a framework applied to spider mites (Acari: Tetranychidae)

V. K. Chaudhari^a,\*^, Abhishek Shukla^b^

^a^ Department of Entomology, N.M. College of Agriculture, Navsari Agricultural University, Navsari 396 450, Gujarat, India

^b^ Department of Entomology, College of Agriculture, Navsari Agricultural University, Waghai Campus 394 730, Gujarat, India

* Corresponding author. E-mail: vkcvaibhav@gmail.com. ORCID: 0000-0001-9688-0134

---

## Highlights

- Reference-library coverage is weighted by whether each label can be verified
- Verifiable coverage of *Tetranychus* climate space is 29% smaller than counts imply
- The gap widens at every coverage threshold tested, from k = 1 to k = 20
- A temporal hold-out shows sequencing lagging the recorded range at six cutoffs
- Workflow takes any taxon name and returns a ranked sampling-priority list

---

## Abstract

Gap analyses of DNA barcode reference libraries count records. Counting assumes that the species label attached to each record is correct, an assumption that fails for taxa where morphological identification is hard and public databases are known to carry misidentifications. We describe a framework that measures instead the *verifiable* coverage of a reference library across the climate space in which a taxon is recorded, and that ranks the localities where one additional reference specimen would add most. Each sequence record carries a weight for whether its label can be confirmed against a morphologically identified anchor specimen; the weight for a record whose cluster holds k anchored species is derived as 1/k rather than chosen. Coverage at a point in climate space is the probability that at least k verifiable sequences lie within one kernel bandwidth of it, and the deficit index is the share of recorded occurrence lying outside that coverage. We applied the framework to spider mites, using 9,915 georeferenced *Tetranychus* occurrence records from GBIF and 974 georeferenced COI records from GenBank and BOLD, with WorldClim 2.1 bioclimatic variables at each point. Weighting reduced 974 sequence records to 405.6 effective sequences and raised the deficit index from 0.218 to 0.282, a relative widening of 29.2%. The widening held between 28.1% and 30.8% across coverage thresholds spanning k = 1 to k = 20, and the weighted index exceeded the unweighted index even when every unverifiable record was credited in full (0.223) and when every weight was replaced by the observed probability that a label is merely correct (0.220). A temporal hold-out trained on records up to a cutoff year placed post-cutoff sequencing below a recorded-range null at all six cutoffs with an adequate training library, with observed-to-null ratios of 0.69 to 0.97: sequencing effort has not been directed at the climatic gaps. *Oligonychus*, for which no anchor verdicts were available, returned an unweighted deficit index of 0.789. The result is not a claim that public labels are wrong — unchecked labels in this dataset stood in 91.2% of checkable cases — but that a library user cannot tell which are right, so unverifiable records cannot carry the weight of a reference. All code, derived data and results are archived and citable.

**Keywords:** DNA barcoding; COI; reference library; coverage gap analysis; Tetranychidae; biodiversity informatics

---

## 1. Introduction

Species identification from short standardised DNA fragments depends on a reference library in which sequences carry correct species names. The size of that dependence has grown with metabarcoding and environmental DNA, where the identity assigned to an unknown sequence is inherited wholesale from its closest reference match. Assessments of whether reference libraries are adequate have accordingly become a recognised genre. Weigand et al. (2019) analysed coverage of 28,000 aquatic species across BOLD and GenBank for European biomonitoring and reported that coverage varies strongly among taxonomic groups and geographic regions, with many groups holding fewer than five barcodes per species. Comparable regional gap analyses have since been published for Atlantic Iberian macroinvertebrates and for annelids of the Arabian Peninsula, the latter also finding that 43% of the relevant BOLD Barcode Index Numbers carried taxonomic ambiguities (Teixeira et al., 2025).

These assessments share a measurement convention: a species is covered when the library holds some number of records bearing its name. The convention treats the name as given. Whether that is safe depends on the taxon. For metazoan mitochondrial sequences as a whole it appears largely safe — Leray et al. (2019) clustered 4.7 million GenBank sequences and estimated an error rate below 1% at the genus level. At lower taxonomic levels and in difficult groups the picture is different. Pentinsaari et al. (2020) traced reported identification failures in BOLD and GenBank largely to errors in sequence acquisition and interpretation rather than to the libraries themselves, while Carné et al. (2025) found errors in GenBank's tetrapod records at the level of the vouchers underpinning them and proposed third-party annotation as a remedy.

Spider mites are a group where the convention is unsafe, and demonstrably so. Morphological separation within *Tetranychus* sensu stricto requires slide-mounted males and a small and declining pool of specialists (Ros and Breeuwer, 2007). *Tetranychus urticae*, *T. cinnabarinus*, *T. turkestani* and *T. truncatus* are not reliably separated by COI (Li et al., 2010), and the species status of *T. cinnabarinus* has since been resolved as a red form of the polymorphic *T. urticae* (Auger et al., 2013). Against that background, de Mendonça et al. (2011) evaluated 138 COI and 105 ITS2 GenBank accessions deposited under *T. urticae* and related taxa and estimated that close to 30% were misidentified or dubious. A coverage figure built by counting records bearing the name *T. urticae* inherits that error rate without recording it.

There is a second, independent problem with counting. Reference sequences are not drawn at random from the places a species occurs. Occurrence records themselves are spatially and environmentally biased, and a library assembled opportunistically from those records compounds the bias. A coverage statement that is taxonomic only — so many species with so many barcodes — says nothing about whether the sequences that exist come from the range of environments in which the species is actually found. For a pest of global distribution, whose populations differ in thermal biology and whose identification matters at ports and in protected cultivation, the environmental dimension of coverage is the one that determines whether the library will work where it is needed.

We set out to measure both at once. The framework described here asks, for a given taxon, how much of the climate space in which it is recorded is backed by reference sequences whose species label can be verified against a morphologically identified specimen, and where a new specimen would most reduce the shortfall. The quantity is deliberately verifiability rather than accuracy. We take no position on how often public labels are right; we measure how much of the library a user could check. Spider mites serve as the demonstration because an independent per-record verification against morphologically identified anchors already exists for them, but every step of the workflow takes a taxon name as an argument and nothing in it is specific to Tetranychidae.

## 2. Materials and methods

### 2.1. Occurrence records

Occurrence records were obtained from GBIF through the download API, which issues a citable DOI for a frozen snapshot, on 31 October 2026: *Tetranychus urticae* (6,646 records; doi:10.15468/dl.xxfcra), *Tetranychus* at genus level (15,793 records; doi:10.15468/dl.aejp4w) and *Oligonychus* at genus level (4,706 records; doi:10.15468/dl.5mckab).

Coordinates were screened before use. Records were removed when they carried a GBIF issue flag that makes a coordinate unusable (`ZERO_COORDINATE`, `COORDINATE_OUT_OF_RANGE`, `COORDINATE_INVALID`, `COUNTRY_COORDINATE_MISMATCH`, `PRESUMED_SWAPPED_COORDINATE`, `PRESUMED_NEGATED_LATITUDE`, `PRESUMED_NEGATED_LONGITUDE`), when latitude equalled longitude, when both were integers, when the stated coordinate uncertainty exceeded 50 km, or when a flag indicated derivation from an administrative centre. Records whose basis was `HUMAN_OBSERVATION` or `MACHINE_OBSERVATION` were separated rather than deleted: a spider mite cannot be identified to species from a photograph, so these were held back from the main analysis and retained for sensitivity testing. For the genus-level *Tetranychus* download this screen retained 9,623 of 15,793 records, removing 3,894 without coordinates, 1,227 flagged as centroid-derived, 133 with fatal issues, 36 with integer coordinates and 26 with excessive uncertainty, and holding back 854 photograph-based observations.

### 2.2. Reference sequence records

GenBank and BOLD records for *Tetranychus* and *Oligonychus* were taken from the archived deposit of a preceding study of these libraries (doi:10.5281/zenodo.23808136), so that the library measured here is the one the identification verdicts describe. Sequence records were georeferenced from the GenBank `lat_lon` qualifier and from BOLD coordinates; BOLD coordinates flagged as centroid-derived were excluded at extraction. Collection year was recovered from the GenBank `collection_date` qualifier, the BOLD collection date fields and the GBIF `year` and `eventDate` fields.

Only the barcode region can supply reference coverage for a COI library. Marker was therefore treated as a filter rather than as an attribute: a georeferenced record whose sequence was COI-3P, another gene, or absent was moved to the occurrence side of the analysis, where it still locates the taxon but no longer counts as reference coverage. For *Tetranychus* this moved 94 COI-3P records, one CYTB record and one 18S record.

### 2.3. Identification verdicts and verifiability weights

Each sequence record carries a per-record verdict from the preceding study, in which public records were tested against morphologically identified anchor specimens drawn from taxonomic publications. The verdicts and their weights are given in Table 1.

**Table 1.** Verifiability weights. Only `AGREES_SHARED_CLUSTER` is derived from the data; the remainder state what an unverifiable record is worth as a reference and are not estimates of how often such labels are correct.

| Verdict | Weight | Meaning |
|---|---|---|
| `AGREES` | 1.00 | label confirmed against an anchored specimen in the same cluster |
| `AGREES_SHARED_CLUSTER` | 1/k | cluster holds k anchored species, so the label fits any one of them |
| `NO_ANCHOR` | 0.25 | no anchored specimen in the cluster; label unverifiable |
| `UNNAMED_SUGGESTED` | 0.25 | record unnamed, one species suggested by clustering |
| `UNNAMED_SEVERAL_POSSIBLE` | 0.10 | record unnamed, several species possible |
| `DISAGREES` | 0.00 | label contradicted by an anchored specimen |

The shared-cluster weight is derived rather than chosen. A record whose cluster contains k anchored species is equally consistent with any of them, so its label is worth 1/k. In the present data every such record sat in a cluster holding exactly two anchored species, giving 0.50 throughout; the rule is stated because k is not a constant in general. Records with no verdict took a default weight of 0.25, and the consequences of that default are bounded explicitly (Section 2.8).

### 2.4. Climate space

WorldClim 2.1 bioclimatic variables BIO1–BIO19 at 10 arc-minute resolution (Fick and Hijmans, 2017) were sampled at the coordinates of every georeferenced record. The 19 variables were standardised and reduced by singular value decomposition to the first two principal components, fitted jointly across occurrence and sequence records so that both sets occupy one space. A 120 × 120 grid was laid over the occupied region with a 5% margin.

### 2.5. Coverage and the deficit index

Occurrence and sequence records were summarised as Gaussian kernel densities on the grid, the sequence density weighted by verifiability. Writing *A* for the effective area of one kernel bandwidth, 2π√det(**H**), and *N* for the sum of verifiability weights, the expected number of verifiable sequences within one bandwidth of a point *z* is

λ(*z*) = *N* · *f*(*z*) · *A*

and coverage saturates as that number approaches a stated threshold *k*:

*C*(*z*) = 1 − exp(−λ(*z*)/*k*)

Because *N* is a sum of weights rather than a count, reducing confidence in the labels lowers coverage in absolute terms rather than redistributing it. The deficit index is the occurrence-density-weighted mean shortfall,

*D* = ∫ *p*(*z*)(1 − *C*(*z*)) d*z*

where *p* is occurrence density normalised to unit mass, so *D* = 0 is full coverage and *D* = 1 is none. We set *k* = 5 for the primary analysis. There is no standard value; five follows the threshold at which (Weigand et al., 2019) reported barcode coverage alongside the single-record threshold, and the index is reported across *k* = 1 to 20.

### 2.6. Sampling priorities

Each occurrence record was assigned the uncovered fraction 1 − *C* at its position in climate space, and records were binned into 2° geographic cells. Cells were ranked by uncovered mass, the sum of uncovered fractions they contain. Selection was greedy with climatic suppression: after a cell was chosen, the uncovered fraction of every remaining record was multiplied by 1 − exp(−½(*d*/*h*)²), where *d* is the distance in principal-component space to the chosen cell's centroid and *h* the kernel bandwidth, so that a specimen from one cell also counts as covering climatically similar cells and the list spreads across climate space rather than selecting neighbours.

### 2.7. Temporal hold-out

To test whether the deficit surface describes something that sampling has responded to, the analysis was trained on records collected up to a cutoff year and tested against sequence records collected afterwards. Two null models were used, both drawing locations from the pre-cutoff occurrence cloud with replacement, 2,000 draws each. The *recorded-range* null weights draws by occurrence density and asks whether sequencing beat sampling the known range at random. The *existing-sequencing* null weights draws by pre-cutoff sequence density and guards the alternative reading that any new locality looks novel against a sparse library whatever its climate. Cutoffs from 2006 to 2020 were swept; those whose training library held fewer than 150 sequences were flagged as too thin for a fair test and are reported but not interpreted.

### 2.8. Sensitivity

Three families of sensitivity analysis were run. Weight schemes spanned a naive scheme in which every record counts fully, strict and lenient variants, an anchored-records-only scheme, a scheme in which every record lacking a verdict is credited in full, and a *label-accuracy floor* in which every weight is replaced by the observed probability that a label is simply correct. That probability was estimated from records whose cluster contained an anchored specimen, as one minus the share of such records contradicted. The coverage threshold *k* was swept from 1 to 20. Each run also reports the index with unverified records credited in full and with them discarded, bounding the effect of the default weight.

### 2.9. Software and data availability

The workflow is implemented in Python 3.12 using NumPy, SciPy, rasterio and Biopython, with 30 automated tests. Code, derived data tables and all results are archived at Zenodo and the repository is public at github.com/vkcvaibhav-eng/reference-coverage. Every script takes a taxon name as an argument; reproducing the analysis for another taxon requires no change to the code.

### 2.10. Use of generative AI

The analysis code, figures and a first draft of this manuscript were produced with the assistance of a large language model (Claude, Anthropic). We designed the study, selected and verified the sources, checked the analysis and its outputs, and revised the text. We are responsible for the content.

## 3. Results

### 3.1. The analysed set

After screening, the *Tetranychus* analysis used 10,889 georeferenced records carrying a complete BIO1–BIO19 vector: 9,915 occurrence records and 974 reference sequence records. The first two principal components of the bioclimatic variables accounted for 47.4% and 24.4% of variance.

Of the 974 sequence records, 688 carried an identification verdict and 286 did not. Mean verifiability was 0.416, so the 974 records amounted to **405.6 effective sequences** — a 2.4-fold reduction between what the library contains and what a user could check. Among all records in the anchor table that could be checked against an anchored specimen, 1,288 of 1,412 labels stood and 124 were contradicted, an observed label accuracy of 0.912.

### 3.2. Coverage deficit

Counting records gave a deficit index of **0.218**. Weighting by verifiability gave **0.282**, a widening of 0.064 or **29.2% in relative terms** (Fig. 1). Under the weighted surface, 745 of 9,915 occurrence records — 7.5% — sat in climate space with less than 10% effective coverage.

The widening was not an artefact of the coverage threshold. Across *k* = 1, 3, 5, 10 and 20 the index rose from 0.086 to 0.420 unweighted and from 0.112 to 0.540 weighted, but the relative widening stayed between 28.1% and 30.8% (Fig. 3b). The absolute level of the index depends on the threshold; the effect of weighting does not.

![Fig1](figs/Fig1_climate_space.png)

**Fig. 1.** Climate space of *Tetranychus*, summarised on the first two principal components of BIO1–BIO19. (a) Occurrence density from 9,915 georeferenced GBIF records. (b) Verifiable coverage from 974 COI records weighted by identification verifiability, at a threshold of five verifiable sequences per kernel bandwidth. (c) The resulting coverage deficit, the product of occurrence mass and the shortfall in coverage. White contours show occurrence density in all three panels.


Nor did it depend on how the 286 records without a verdict were treated. Crediting them in full gave 0.223 and discarding them gave 0.305, both above the unweighted 0.218. The label-accuracy floor — every weight replaced by the observed 0.912 chance that a label is merely correct — gave 0.220, still above 0.218 but only marginally (Fig. 3a). An anchored-records-only scheme, retaining 132 effective sequences, gave 0.378.

![Fig3](figs/Fig3_sensitivity.png)

**Fig. 3.** Sensitivity of the deficit index. (a) Weight schemes, with effective sequence counts in parentheses; the dashed line is the unweighted index of 0.218. (b) Coverage threshold *k* swept from 1 to 20, with the relative widening between the counted and weighted indices annotated at each value.


### 3.3. Where the gaps are

The 9,915 occurrence records fell in 655 two-degree cells. Ranked by uncovered mass, the leading cells were in northern Australia (12–16°S, 130–132°E; 305 occurrence records, mean uncovered fraction 0.640), south-eastern Australia (32–36°S, 150–152°E; 1,191 records, 0.147), Réunion (20–24°S, 54–56°E; 195 records, 0.593), Seychelles (4–8°S, 54–56°E; 75 records, 0.863), central China (36–40°N, 104–106°E; 53 records, 0.735) and central Japan (34–38°N, 136–138°E; 39 records, 0.623) (Fig. 2).

![Fig2](figs/Fig2_geographic_deficit.png)

**Fig. 2.** Coverage deficit mapped back to geography. Each point is a 2° cell holding *Tetranychus* occurrence records; colour gives the mean uncovered fraction of the records it contains and symbol size scales with record count. Open circles mark the ten highest-ranked sampling priorities by uncovered mass, numbered for the first five.


The two ranking criteria separate two kinds of gap. Uncovered mass favours cells holding many occurrence records, so south-eastern Australia ranks second on 1,191 records despite a mean uncovered fraction of only 0.147. Mean uncovered fraction favours cells whose climate is unrepresented whatever their record count, which is why Seychelles, on 75 records, carries the highest value in the top six at 0.863.

### 3.4. Temporal hold-out

At the 2015 cutoff, the 519 training sequences gave a deficit index of 0.523, and the 354 sequence records collected afterwards fell at a mean uncovered fraction of 0.366. The recorded-range null expected 0.414 (95% interval 0.373–0.452; one-sided *p* = 0.991) and the existing-sequencing null expected 0.095 (0.079–0.112; *p* = 0.0005).

The pattern was consistent across the sweep (Fig. 4). At all six cutoffs with a training library of at least 150 sequences — 2012, 2014, 2015, 2016, 2018 and 2020 — observed uncovered fraction fell **below** the recorded-range null, with ratios of 0.811, 0.844, 0.884, 0.970, 0.859 and 0.694. At the three thin cutoffs (2006, 2008, 2010; 34 to 90 training sequences) the ratio exceeded one, but those training libraries are too small to constitute a fair test. Against the existing-sequencing null the observed value exceeded expectation at every cutoff (*p* = 0.0005 throughout).


![Fig4](figs/Fig4_holdout.png)

**Fig. 4.** Temporal hold-out. (a) Mean uncovered fraction of post-cutoff sequence records against both nulls, by cutoff year; the shaded band marks cutoffs whose training library held fewer than 150 sequences. (b) The same result as a ratio to the recorded-range null, with the training library size printed below each point. Values below one indicate sequencing falling in less uncovered climate than draws from the recorded range.

Post-cutoff sequencing therefore moved away from where sequencing had already been, but not towards the climates the library lacked: relative to the recorded range it lagged.

### 3.5. *Oligonychus*

The *Oligonychus* analysis used 3,199 records (2,878 occurrence, 321 sequence). **No sequence record carried an identification verdict**, the anchor table of the preceding study covering none of the georeferenced *Oligonychus* records. Its weighting is therefore entirely the default value and carries no verification information. We report the unweighted index only: **0.789**, against 0.218 for *Tetranychus*. At the 2015 cutoff the training library held 23 sequences, far below the threshold for a fair hold-out, so no temporal test is reported.

## 4. Discussion

### 4.1. What the weighting changes, and what it does not

The central result is narrow and should be stated narrowly. Weighting by verifiability widened the measured gap in *Tetranychus* COI coverage by about 29%, and that widening was stable across every coverage threshold and every treatment of the unverified records that we tested. It is not a finding that public labels are wrong. In this dataset 91.2% of checkable labels stood, close to the general picture for metazoan mitochondrial sequences (Leray et al., 2019) and consistent with the view that reported identification failures often originate in sequence handling rather than in the libraries (Pentinsaari et al., 2020).

The label-accuracy floor makes the boundary explicit. When every weight is replaced by the bare probability that a label is correct, the index falls to 0.220 against an unweighted 0.218 — the effect all but disappears. A framework that claimed to measure accuracy would collapse at that point. What survives is the claim we make: that a library user cannot distinguish a correct unanchored label from an incorrect one, so unverifiable records cannot carry the weight of a reference. The distinction matters for how the result should be read and for how it should be cited.

One caution attaches to the floor itself. The 0.912 base rate was estimated only from records whose cluster contained an anchored specimen, and those clusters are the better-studied part of the library. Transferring that rate to records with no anchor is an optimistic extrapolation rather than a neutral one, which is why we report it as a floor and not as a working weight.

### 4.2. Verifiability as a coverage criterion

Existing gap analyses count species with at least one, or at least five, barcode records (Weigand et al., 2019). The threshold convention is sound and we have adopted it, but it is applied to records whose labels are taken as given. For taxa where identification is contested the convention produces a coverage figure that inherits an unrecorded error rate. In *Tetranychus* the rate is not small: de Mendonça et al. (2011) put misidentified or dubious sequences at close to 30% of the accessions they examined, in a genus where the barcode region does not separate *T. urticae* from several congeners (Li et al., 2010) and where one nominal species has since been synonymised with another (Auger et al., 2013). Counting records in that setting measures the size of the library, not its usable extent.

Making verifiability an explicit weight has a practical consequence beyond the headline number. Because the weight is derived from an independent per-record check against morphologically identified specimens, it carries the taxonomic effort forward into the coverage statistic rather than leaving it in a separate quality assessment. Assessments that quantify library ambiguity and assessments that quantify library coverage have so far run on parallel tracks — Teixeira et al. (2025), for instance, report both a coverage gap and a 43% rate of BIN-level taxonomic ambiguity, but as separate findings. The framework here combines them into one quantity.

### 4.3. Climate space as the unit of coverage

Measuring coverage over climate space rather than over a species checklist changes what a gap means. A checklist gap says a species has no barcode. A climate-space gap says that populations living under conditions represented in the occurrence record have no verifiable reference sequence, even where the species as a whole is well sequenced. For *Tetranychus*, which is sequenced in the hundreds, the checklist question is already answered; the climate-space question is not.

The cells our ranking returns are interpretable in those terms. Réunion and Seychelles carry high mean uncovered fractions on modest record counts, consistent with island climates under-represented in a library assembled largely in continental temperate agriculture. Northern Australia combines a high uncovered fraction with a substantial record count, which is why it leads on uncovered mass. South-eastern Australia ranks high on mass alone: its climate is comparatively well covered, and it earns its place through the sheer number of occurrence records it holds.

That last case exposes a limitation. Uncovered mass inherits whatever bias the occurrence record carries, and GBIF occurrence data are not an unbiased sample of where a species lives. Our ranking is a statement about two public evidence layers, the occurrence record and the reference library, and about the relation between them. It is not a statement about the true distribution of the mites. Readers wanting a ranking less sensitive to collection effort should use the mean uncovered fraction column, which is reported alongside uncovered mass for every cell.

### 4.4. The hold-out result

We designed the temporal hold-out expecting it to validate the deficit surface by showing that subsequent sequencing had been drawn towards the gaps. It did not. At every cutoff with an adequate training library, post-cutoff sequence records fell in less uncovered climate than draws from the recorded range, with ratios between 0.69 and 0.97.

Two readings are available and the data do not separate them. Sequencing may be following the taxon's economic importance, or accessible collections, or existing collaborations, none of which track climatic novelty. Alternatively the deficit surface may be identifying climates that are genuinely hard to sample — remote, politically difficult, or supporting populations that are rarely collected for reasons unrelated to their interest. The contrast with the existing-sequencing null is informative: post-cutoff records did fall in more uncovered climate than draws weighted by where sequencing had already been, at every cutoff. New sequencing has moved, but it has moved more slowly than the recorded range has expanded.

Either way, the hold-out is no longer a validation that the metric predicts sampling behaviour, and we do not present it as one. It establishes something weaker and still useful: the deficit surface is not merely a map of where people work. Had it been, new sequencing would have landed in high-deficit cells automatically.

### 4.5. *Oligonychus* and the limits of the framework

The *Oligonychus* result is a demonstration of a precondition rather than a finding about the genus. Its unweighted deficit index of 0.789 is high, but with no verdict attached to any of its 321 georeferenced sequence records we cannot say what share of that coverage is verifiable. The framework requires an independent per-record verification as input. Where none exists it reduces to a conventional, count-based gap analysis over climate space — useful, but not what we set out to measure. Extending anchor-based verification to *Oligonychus* is the obvious next step and would make the two genera comparable.

### 4.6. Other limitations

WorldClim 2.1 is a 1970–2000 climatology (Fick and Hijmans, 2017), and the records analysed span decades; a record collected in 1960 and one collected in 2020 receive the same climate values. At 10 arc-minute resolution a protected-cultivation site and open field 15 km away are indistinguishable, which matters for a pest much of whose distribution lies under glass and plastic. We have not separated protected from open cultivation, and a land-cover or plastic-covered-agriculture layer would allow it.

We have also not tested whether the gaps track accessibility and research capacity rather than climate. Occurrence records cluster near roads and institutions, and a deficit surface built on them will inherit that structure. Modelling the deficit against travel time and per-country research output alongside climate would separate the explanations; we regard this as the most important extension of the present work.

Finally, the default weight of 0.25 applied to the 286 *Tetranychus* records without a verdict is a stated position, not a measurement. We have bounded its influence rather than removed it. Extending the verification to those records would remove the assumption entirely.

## 5. Conclusions

Public COI reference libraries are routinely assessed by counting records, which presumes that each record's species label is correct. For spider mites that presumption is known to fail, and the framework presented here replaces it with an explicit weight for whether a label can be verified against a morphologically identified specimen. Applied to *Tetranychus*, 974 georeferenced COI records resolved to 405.6 effective sequences and the climate-space coverage deficit rose from 0.218 to 0.282, a 29.2% relative widening that held across coverage thresholds from k = 1 to k = 20 and survived both bounds on the treatment of unverified records. A temporal hold-out showed post-cutoff sequencing falling below a recorded-range null at all six cutoffs with an adequate training library, indicating that sequencing effort has lagged rather than closed the climatic gaps. The measured quantity is verifiable coverage, not accuracy: unchecked labels in this dataset stood in 91.2% of checkable cases, but a library user cannot tell which, and that is precisely why the distinction is worth making. Every script takes a taxon name as an argument, and the framework applies to any group for which a per-record verification against morphologically identified specimens exists.

---

## Declarations

**Funding.** No funds, grants or other support were received for conducting this study.

**Competing interests.** The authors declare no competing interests.

**CRediT author statement.** **V. K. Chaudhari:** Conceptualization, Methodology, Software, Formal analysis, Data curation, Writing – original draft, Writing – review & editing, Visualization. **Abhishek Shukla:** Conceptualization, Validation, Writing – review & editing, Supervision.

**Data availability.** Code, derived data and results are archived at Zenodo and available at github.com/vkcvaibhav-eng/reference-coverage. Occurrence data are the GBIF downloads doi:10.15468/dl.xxfcra, doi:10.15468/dl.aejp4w and doi:10.15468/dl.5mckab. Identification verdicts are from doi:10.5281/zenodo.23808136.

---

## References

Auger, P., Migeon, A., Ueckermann, E.A., Tiedt, L., Navajas, M., 2013. Evidence for synonymy between *Tetranychus urticae* and *Tetranychus cinnabarinus* (Acari, Prostigmata, Tetranychidae): review and new data. Acarologia 53, 383–415.

Carné, A., Vieites, D.R., van den Burg, M.P., 2025. In vouchers we (hope to) trust: unveiling hidden errors in GenBank's tetrapod taxonomic foundations. Molecular Ecology 34, e17812. https://doi.org/10.1111/mec.17812

de Mendonça, R.S., Navia, D., Diniz, I.R., Auger, P., Navajas, M., 2011. A critical review on some closely related species of *Tetranychus* sensu stricto (Acari: Tetranychidae) in the public DNA sequences databases. Experimental and Applied Acarology 55, 1–23. https://doi.org/10.1007/s10493-011-9453-5

Fick, S.E., Hijmans, R.J., 2017. WorldClim 2: new 1-km spatial resolution climate surfaces for global land areas. International Journal of Climatology 37, 4302–4315. https://doi.org/10.1002/joc.5086

Leray, M., Knowlton, N., Ho, S.-L., Nguyen, B.N., Machida, R.J., 2019. GenBank is a reliable resource for 21st century biodiversity research. Proceedings of the National Academy of Sciences 116, 22651–22656. https://doi.org/10.1073/pnas.1911714116

Li, G.-Q., Xue, X.-F., Zhang, K.-J., Hong, X.-Y., 2010. Identification and molecular phylogeny of agriculturally important spider mites (Acari: Tetranychidae) based on mitochondrial and nuclear ribosomal DNA sequences, with an emphasis on *Tetranychus*. Zootaxa 2647, 1–15.

Pentinsaari, M., Ratnasingham, S., Miller, S.E., Hebert, P.D.N., 2020. BOLD and GenBank revisited — do identification errors arise in the lab or in the sequence libraries? PLoS ONE 15, e0231814. https://doi.org/10.1371/journal.pone.0231814

Ros, V.I.D., Breeuwer, J.A.J., 2007. Spider mite (Acari: Tetranychidae) mitochondrial COI phylogeny reviewed: host plant relationships, phylogeography, reproductive parasites and barcoding. Experimental and Applied Acarology 42, 239–262. https://doi.org/10.1007/s10493-007-9092-z

Teixeira, M.A.L., et al., 2025. Gaps and data ambiguities in DNA reference libraries: a limiting factor for molecular-based biodiversity assessments using annelids as a case study. Ecology and Evolution 15, e71544. https://doi.org/10.1002/ece3.71544

Weigand, H., Beermann, A.J., Čiampor, F., Costa, F.O., Csabai, Z., Duarte, S., Geiger, M.F., Grabowski, M., Rimet, F., Rulik, B., Strand, M., Szucsich, N., Weigand, A.M., Willassen, E., Wyler, S.A., Bouchez, A., Borja, A., Čiamporová-Zaťovičová, Z., Ferreira, S., Dijkstra, K.-D.B., Eisendle, U., Freyhof, J., Gadawski, P., Graf, W., Haegerbaeumer, A., van der Hoorn, B.B., Japoshvili, B., Keresztes, L., Keskin, E., Leese, F., Macher, J.N., Mamos, T., Paz, G., Pešić, V., Pfannkuchen, D.M., Pfannkuchen, M.A., Price, B.W., Rinkevich, B., Teixeira, M.A.L., Várbíró, G., Ekrem, T., 2019. DNA barcode reference libraries for the monitoring of aquatic biota in Europe: gap-analysis and recommendations for future work. Science of the Total Environment 678, 499–524. https://doi.org/10.1016/j.scitotenv.2019.04.247
