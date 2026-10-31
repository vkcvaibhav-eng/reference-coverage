# Evidence-weighted climatic coverage of public COI reference libraries: a framework applied to spider mites

V. K. Chaudhari^a,\*^, Abhishek Shukla^b^

^a^ Department of Entomology, N.M. College of Agriculture, Navsari Agricultural University, Navsari 396 450, Gujarat, India

^b^ Department of Entomology, College of Agriculture, Navsari Agricultural University, Waghai Campus 394 730, Gujarat, India

\* Corresponding author. E-mail: vkcvaibhav@gmail.com. ORCID: 0000-0001-9688-0134

---

## Highlights

- Reference coverage is weighted by whether each species label can be verified
- Duplicate identifiers accounted for 25.3% of Tetranychus database rows
- Verifiability raises the Tetranychus coverage deficit from 0.237 to 0.305
- Candidate localities are ranked by one specimen's exact effect on that index
- Later-collected records occupy better-covered climate than the recorded range

---

## Abstract

Gap analyses of DNA barcode reference libraries can overstate coverage when database identifiers are duplicated and the taxonomic evidence behind a label is incomplete. We measure instead the *verifiable* coverage of a reference library across the climate space in which a taxon is recorded. Records are first resolved to specimens, since GenBank and BOLD index the same material under different identifiers. Each specimen carries a weight for whether its label can be confirmed against a morphologically identified anchor, with an equal-ambiguity convention of 1/k for a cluster holding k anchored species. The deficit index is the occurrence-density-weighted mean shortfall in the coverage score; localities are ranked by the exact fall in that index when one verified specimen is added. For *Tetranychus*, 974 georeferenced COI rows resolved to 728 specimens, 25.3% of rows being redundant, and resolution raised the verdict-match rate from 70.6% to 94.8%. Weighting left 344.6 effective specimens and raised the deficit index from 0.237 to 0.305; the weighted index exceeded the counted index at every coverage scale from k = 1 to k = 20, by 28.6% to 35.1%. *Oligonychus*, 44.6% of whose rows were redundant, gave 0.703 and 0.745. A retrospective split on collection year placed later-collected records in significantly better-covered climate than the recorded range. The quantity is verifiability, not accuracy: 90.9% of checkable labels were not contradicted, but most sit in clusters that do not resolve to one species.

**Keywords:** DNA barcoding; COI; reference library; coverage; Tetranychidae; biodiversity informatics

---

## 1. Introduction

Species identification from short standardised DNA fragments depends on a reference library in which sequences carry correct species names. Metabarcoding and environmental DNA have deepened that dependence, because the identity assigned to an unknown sequence is inherited from its closest reference match. Assessments of whether reference libraries are adequate have accordingly become a recognised genre. Weigand et al. (2019) analysed coverage of 28,000 aquatic species across BOLD and GenBank for European biomonitoring and reported that coverage varies strongly among taxonomic groups and geographic regions, with many groups holding fewer than five barcodes per species. Comparable regional analyses have since been published for annelids of the Arabian Peninsula, which also found that 43% of the relevant Barcode Index Numbers carried taxonomic ambiguities (Teixeira et al., 2025).

These assessments share a measurement convention: a species is covered when the library holds some number of records bearing its name. The convention treats both the record and the name as given. Whether that is safe depends on the taxon. For metazoan mitochondrial sequences as a whole the names appear largely safe — Leray et al. (2019) clustered 4.7 million GenBank sequences and estimated an error rate below 1% at the genus level. In difficult groups, and at lower taxonomic levels, the picture differs. Pentinsaari et al. (2020) traced reported identification failures in BOLD and GenBank largely to errors in sequence acquisition and interpretation rather than to the libraries themselves, while Carné et al. (2025) found errors in GenBank's tetrapod records at the level of the vouchers underpinning them.

Spider mites are a group where the convention is demonstrably unsafe. Morphological separation within *Tetranychus* sensu stricto requires slide-mounted males and a small and declining pool of specialists (Ros and Breeuwer, 2007). *Tetranychus urticae*, *T. cinnabarinus*, *T. turkestani* and *T. truncatus* are not reliably separated by COI (Li et al., 2010), and the species status of *T. cinnabarinus* has since been resolved as a red form of the polymorphic *T. urticae* (Auger et al., 2013). Against that background, de Mendonça et al. (2011) evaluated 138 COI and 105 ITS2 GenBank accessions deposited under *T. urticae* and related taxa and judged close to 30% misidentified or dubious.

The record itself is also less solid than it appears. BOLD deposits its sequences in GenBank, so one mite commonly carries two identifiers and appears in a pooled analysis twice. We found this to be not a marginal effect but a quarter of the *Tetranychus* rows and nearly half the *Oligonychus* rows, and it cuts both ways: duplicate rows inflate the apparent size of a library, while an identifier mismatch can make existing taxonomic evidence appear to be missing.

A second, independent problem attends counting. Reference sequences are not drawn at random from the places a taxon occurs, and occurrence records are themselves spatially and environmentally biased. A coverage statement that is taxonomic only says nothing about whether the sequences that exist come from the range of environments in which the taxon is found. For pests of global distribution, whose populations differ in thermal biology and whose identification matters at ports and in protected cultivation, the environmental dimension determines whether the library will work where it is needed.

Quantifying occurrence in environmental rather than geographic space is established practice. Broennimann et al. (2012) set out a framework applying kernel smoothers to occurrence densities in gridded environmental space, deliberately so that the result does not depend on sampling effort or on an arbitrary resolution, and used it to compare niches between taxa. We adopt that geometry. What we add is not the environmental-space representation but the evidential weighting within it, and the resolution of records to specimens that must precede it.

The framework therefore asks, for a given taxon, how much of its recorded climate space is covered by reference sequences whose label can be verified against a morphologically identified specimen, and where one further specimen would most reduce the shortfall. Coverage is weighted by that evidence rather than counted. The quantity is deliberately verifiability rather than accuracy: we take no position on how often public labels are right, and measure how much of the library a user could check.

## 2. Materials and methods

### 2.1. Occurrence records

Occurrence records were obtained from GBIF through the download API, which issues a citable DOI for a frozen snapshot, on 31 October 2026: *Tetranychus urticae* (6,646 records; doi:10.15468/dl.xxfcra), *Tetranychus* at genus level (15,793 records; doi:10.15468/dl.aejp4w) and *Oligonychus* at genus level (4,706 records; doi:10.15468/dl.5mckab).

Coordinates were screened following the categories of problem catalogued by Zizka et al. (2019). Records were removed when they carried a GBIF issue flag that makes a coordinate unusable (`ZERO_COORDINATE`, `COORDINATE_OUT_OF_RANGE`, `COORDINATE_INVALID`, `COUNTRY_COORDINATE_MISMATCH`, `PRESUMED_SWAPPED_COORDINATE`, `PRESUMED_NEGATED_LATITUDE`, `PRESUMED_NEGATED_LONGITUDE`), when latitude equalled longitude, when both were integers, when stated coordinate uncertainty exceeded 50 km, or when a flag indicated derivation from an administrative centre. These are operational filters on the fields GBIF supplies, not determinations about how a record was made; `COORDINATE_ROUNDED` in particular indicates decimal rounding rather than an established administrative centroid, and we treat it conservatively.

Records whose basis was `HUMAN_OBSERVATION` or `MACHINE_OBSERVATION` were excluded from the main analysis and written to separate files that the extraction step does not read. The label denotes an observation made by a person, not that the record is photograph-based, and an expert field determination may legitimately carry it; we exclude these records because the basis field alone does not establish that a specimen was examined. For the genus-level *Tetranychus* download the screen retained 9,623 of 15,793 records, removing 3,894 without coordinates, 1,227 with centroid-indicating flags, 133 with fatal issues, 36 with integer coordinates and 26 with excessive uncertainty, and setting aside 854 observations.

### 2.2. Reference records and resolution to specimens

GenBank and BOLD records for *Tetranychus* and *Oligonychus* were taken from the archived deposit of a preceding study of these libraries (doi:10.5281/zenodo.23808136), so that the library measured here is the one the identification verdicts describe. Records were georeferenced from the GenBank `lat_lon` qualifier and from BOLD coordinates; BOLD coordinates flagged as centroid-derived were excluded at extraction. Collection year was recovered from the GenBank `collection_date` qualifier, the BOLD collection date fields and the GBIF `year` and `eventDate` fields.

A record is not a specimen. We built an alias map from the per-genus crosswalks published by the preceding study, which carry accession and process ID side by side, and resolved every row to a canonical specimen through 4,502 aliases. For *Tetranychus* this collapsed 246 of 974 rows, 25.3% of the total, leaving 728 specimens; for *Oligonychus* it collapsed 143 of 321 rows, 44.6%, leaving 178. Among rows sharing a canonical key, the retained row is the one carrying a verdict, then by source preference, then by identifier, so the outcome does not depend on input order. We refer to the resolved units as specimens throughout; strictly they are canonical records, and establishing that each corresponds to one physical mite would require per-sequence voucher provenance that the deposited metadata does not supply.

Two representations of one specimen disagreed on coordinates in 101 *Tetranychus* and 100 *Oligonychus* groups. At 10 arc-minute resolution only four of the *Tetranychus* groups yield a different BIO1–BIO19 vector; the remainder are differences too small to change the climate inputs. All conflict groups are listed in the released results.

Only the barcode region can supply reference coverage for a COI library. A georeferenced record whose sequence was COI-3P, another gene, or absent was moved to the occurrence side, where it still locates the taxon but does not count as reference coverage: 94 COI-3P, one CYTB and one 18S record in *Tetranychus*. Marker assignment rests on the deposited gene or marker field. For records naming COI without specifying a region, position within COI was not confirmed by alignment, and the accepted set therefore includes sequences whose overlap with the barcode region is unestablished. This is a limitation on eligibility, not a verified property of the retained set.

### 2.3. Verdicts and verifiability weights

Each specimen carries a per-record verdict from the preceding study, in which public records were tested against morphologically identified anchor specimens drawn from taxonomic publications. Verdicts were matched through the alias map, which raised the match rate for *Tetranychus* from 70.6% of rows to 94.8% of specimens. Table 1 gives the weights and the resulting composition.

**Table 1.** Verifiability weights and the resulting composition of the *Tetranychus* reference set. Generated from the analysis output.

| Verdict | Weight | Meaning | *Tetranychus* specimens |
|---|---|---|---|
| Agrees with anchor | 1.00 | label concordant with an anchored specimen in the same cluster | 132 |
| Agrees, shared cluster | 1/2 = 0.50 | cluster holds k anchored species, so the label fits any one | 281 |
| No anchor | 0.25 | no anchored specimen in the cluster; label unverifiable | 235 |
| Unnamed, one suggested | 0.25 | record unnamed, one species suggested by clustering | 2 |
| Unnamed, several possible | 0.10 | record unnamed, several species possible | 33 |
| Contradicted | 0.00 | label contradicted by an anchored specimen | 7 |
| No verdict | 0.25 | not covered by the anchor check | 38 |
| **Total** | | | **728** |

A distinction matters in reading Table 1. Only 38 specimens lack a verdict altogether; a further 235 have one, `No anchor`, which records that no anchored specimen was available in their cluster. The first is missing information about the record, the second is information that the record is unverifiable.

The shared-cluster weight is an equal-ambiguity convention, not a calibrated probability. A specimen whose cluster contains k anchored species is consistent with any of them, so we credit 1/k; cluster membership supplies the candidate count, but equal evidential support for each candidate is a modelling choice. In the present data every such specimen sat in a cluster holding exactly two anchored species. *Oligonychus* has no specimens in this category.

Concordance with an anchor is also weaker evidence than re-examination of the specimen itself. It establishes that the label is consistent with the taxonomic determination of another specimen in the same sequence cluster, under the clustering threshold, alignment region and anchor-eligibility rules of the preceding study.

### 2.4. Climate space

WorldClim 2.1 bioclimatic variables BIO1–BIO19 at 10 arc-minute resolution (Fick and Hijmans, 2017) were sampled at the coordinates of every georeferenced record. The 19 variables were standardised and reduced by singular value decomposition to the first two principal components, fitted jointly across occurrence and reference records so that both occupy one space; for *Tetranychus* these accounted for 47.7% and 25.0% of variance, leaving 27.3% outside the plane in which coverage is evaluated. A 120 × 120 grid was laid over the occupied region with a 5% margin.

### 2.5. Coverage and the deficit index

Occurrence and reference records were summarised as Gaussian densities on the grid, the reference density weighted by verifiability. Smoothing covariance is held fixed across weight schemes. This matters: a kernel density estimator that sets its bandwidth from the weighted effective sample size spreads credit differently as well as awarding less of it, so weighted and unweighted surfaces would differ in smoothing as well as in evidence. The covariance is computed once from the unweighted reference points by Scott's rule and reused for every scheme.

Writing *A* for the area of one kernel bandwidth, 2π√det(**H**), and *N* for the sum of verifiability weights, the local kernel-weighted reference support at climate position *z* is

λ(*z*) = *N* · *f*(*z*) · *A*

and the coverage score saturates as that number approaches a stated scale *k*:

*C*(*z*) = 1 − exp(−λ(*z*)/*k*)

Because *N* is a sum of weights rather than a count, reducing confidence in labels lowers the score in absolute terms rather than redistributing it. Two qualifications attach. λ is a sum of distance-decaying kernel similarities, not a count of specimens inside a hard boundary, so it is a measure of local support rather than a tally; and at λ = *k* the score is 0.632, so *k* is a saturation scale rather than a threshold. The deficit index is the occurrence-weighted mean shortfall,

*D* = ∫ *p*(*z*)(1 − *C*(*z*)) d*z*

where *p* is occurrence density normalised to unit mass, so *D* = 0 is full coverage and *D* = 1 is none. *D* is a smooth occurrence-weighted shortfall score, not a measured fraction of populations lacking a barcode. We set *k* = 5 for the primary analysis, following the threshold at which Weigand et al. (2019) report barcode coverage alongside the single-record threshold, and report *D* across *k* = 1 to 20.

Coverage is computed within a genus. A *T. ludeni* reference therefore contributes to the score at climates where *T. urticae* occurs. The index measures the climatic representation of a genus-wide verified collection against a genus-wide occurrence record; it does not establish that any named species has species-specific reference support, and we draw no species-level conclusions from it.

### 2.6. Sampling priorities

Occurrence records were binned into 2° cells by the floor of their coordinates, so a cell labelled −24, 166 spans 24° to 22° S and 166° to 168° E. For each cell we compute the exact fall in *D* produced by adding one verified specimen of unit weight at that cell's climate centroid. Because one added specimen raises λ by exactly *A* times a single kernel evaluated at the centroid, the new index is available in closed form. The quantity is the single-addition gain against the current library on the fixed grid, principal-component representation and kernel covariance. It ranks alternative single additions; it is not an optimal collecting plan for the whole list, which would require recomputing the gain after each selection. Uncovered mass and mean uncovered fraction are reported alongside it for every cell.

### 2.7. Retrospective comparison by collection year

To test whether the deficit surface describes something later collecting has tracked, the analysis was trained on specimens collected up to a cutoff year and compared against those collected afterwards. The split is on collection date, not deposition or publication date: a mite collected in 2012 and sequenced in 2024 falls in the training set. The comparison is retrospective with respect to collecting and cannot speak to how sequencing effort changed over time. The principal-component representation and the anchor information are fitted on the full dated set, so the comparison is a retrospective description rather than a prospective prediction.

Two nulls were used, each resampling with equal probability, 2,000 draws, the uncovered fractions of a stated reference set: the pre-cutoff occurrence records (*recorded range*) and the pre-cutoff reference specimens (*already sequenced*). Weighting such draws by a density estimated from the same points would sample the square of that density and would correspond to neither benchmark. Both tails are reported with the Monte Carlo correction (*r* + 1)/(*m* + 1). Cutoffs from 2006 to 2020 were swept; those whose training library held fewer than 150 specimens are flagged and not interpreted. That threshold is a stated adequacy convention, not a derived requirement. Cutoffs share records, so the sweep describes stability rather than providing independent replication, and records from one study or locality are not resampled as blocks.

### 2.8. Sensitivity

Weight schemes spanned a counted scheme in which every specimen contributes fully, strict and lenient variants, an anchored-only scheme, a scheme crediting records without matched verdicts in full, and a uniform non-contradiction scenario in which every weight is replaced by the observed rate at which checkable labels were not contradicted. That last is a sensitivity scenario, not a demonstrated lower bound. The coverage scale *k* was swept from 1 to 20. Each run also reports the index with those records credited in full and discarded, which bounds the treatment of those 38 reference records without matched verdicts, not all uncertainty in the weights.

### 2.9. Software and data availability

The workflow is implemented in Python 3.12 using NumPy, SciPy, rasterio and Biopython, with 30 automated tests. The repository provides the analysis code, automated tests, figure-generation scripts and machine-readable numerical summaries; Table 1 and the figures are emitted by those scripts, and the numerical statements in the text were cross-checked against the summaries. `fetch_gbif.py`, `deficit.py` and `holdout.py` take a taxon name as an argument; the remaining scripts operate on whatever the preceding steps produced. Applying the framework to another taxon requires GBIF and sequence inputs for that taxon and an equivalent per-record verification, not changes to the code.

## 3. Results

### 3.1. The analysed set

After screening and resolution to specimens, the *Tetranychus* analysis used 9,789 georeferenced records carrying a complete BIO1–BIO19 vector: 9,061 on the occurrence side, comprising 8,965 GBIF records and 96 reclassified non-barcode reference records, and 728 reference specimens. The first two principal components accounted for 47.7% and 25.0% of bioclimatic variance.

Of the 728 specimens, 690 carried a verdict and 38 did not. Mean verifiability was 0.473, so 728 specimens amounted to **344.6 effective specimens**. Across all records in the anchor tables that could be checked against an anchored specimen, 1,563 of 1,720 labels were not contradicted and 157 were, a non-contradiction rate of 0.909. This is not a measured species-label accuracy: 829 of the uncontradicted labels sit in clusters holding more than one anchored species and are therefore consistent with, but not resolved to, the deposited name.

### 3.2. Coverage deficit

Counting specimens gave a deficit index of **0.237**. Weighting by verifiability gave **0.305**, a widening of 0.068 or **28.6% in relative terms** (Fig. 1). Under the weighted surface, 819 of 9,061 occurrence records — 9.0% — sat in climate space with a coverage score below 0.10.
![Fig1](figs/Fig1_climate_space.png)

**Fig. 1.** Climate space of *Tetranychus* on the first two principal components of BIO1–BIO19. (a) Occurrence density from 9,061 georeferenced records, of which 8,965 are GBIF occurrences and 96 are reference records reclassified as non-barcode. (b) Verifiable coverage score from 728 specimens weighted by verifiability, at a saturation scale of five verifiable specimens per kernel bandwidth. (c) The resulting coverage deficit, the product of occurrence mass and the shortfall in the coverage score. White contours show occurrence density in all three panels.


The weighted deficit exceeded the counted deficit at every tested coverage scale, with relative increases ranging from 28.6% to 35.1%; across *k* = 1, 3, 5, 10 and 20 the index itself rose from 0.095 to 0.447 counted and 0.128 to 0.576 weighted (Fig. 2b).

The 38 records without a matched verdict make little difference. Crediting them in full gave 0.284 and discarding them gave 0.315, both above the counted 0.237. The uniform non-contradiction scenario, in which every weight becomes 0.909, gave 0.249 (Fig. 2a). An anchored-only scheme, retaining 132 effective specimens, gave 0.476.
![Fig2](figs/Fig3_sensitivity.png)

**Fig. 2.** Sensitivity of the deficit index. (a) Weight schemes, with effective specimen counts in parentheses; the dashed line is the counted index of 0.237. (b) The coverage saturation scale *k* swept from 1 to 20, with the relative widening between the counted and weighted indices annotated at each value. Relative increases range from 28.6% to 35.1% while the index itself more than quadruples.


### 3.3. Where the gaps are

The 9,061 occurrence records fell in 611 two-degree cells. Ranked by the exact reduction in *D* from one verified specimen, the leading cells were New Caledonia (24–22° S, 166–168° E; Δ*D* = 0.00636 on 136 occurrence records), western Japan (32–34° N, 132–134° E; 0.00627 on 6 records), Réunion (22–20° S, 54–56° E; 0.00623 on 177 records), the United States Gulf coast (30–32° N, 92–90° W; 0.00622 on 8 records) and south-eastern Australia (34–32° S, 152–154° E; 0.00621 on 51 records) (Fig. 3).
![Fig3](figs/Fig2_geographic_deficit.png)

**Fig. 3.** Coverage deficit mapped back to geography. Each point is a 2° cell holding *Tetranychus* occurrence records, plotted at its centre; colour gives the mean uncovered fraction of the records it contains and symbol size scales with record count. Open circles mark the ten cells where one verified specimen would most reduce the deficit index, numbered for the first five. This ranking is not a ranking by uncovered fraction.


The ranking is not a ranking by shortfall. The second cell holds six occurrence records and an uncovered mass of 2.37; the third holds 177 records and a mass of 158.6. What separates them is where their climate sits relative to the existing verified collection, not how much uncovered occurrence they contain, and the gains among the leading cells are close together.

### 3.4. Retrospective comparison by collection year

At the 2015 cutoff, 354 training specimens gave a deficit index of 0.565, and the 289 specimens collected afterwards fell at a mean uncovered fraction of 0.420. A resample of the pre-cutoff occurrence records expected 0.546 and a resample of the pre-cutoff reference specimens expected 0.109.

Across the five cutoffs with a training library of at least 150 specimens — 2014, 2015, 2016, 2018 and 2020 — later-collected records fell in **significantly less** uncovered climate than a resample of the recorded range, with lower-tail *p* = 0.0005 at every cutoff and ratios to that null of 0.775, 0.769, 0.806, 0.674 and 0.540 (Fig. 4). Against the already-sequenced null the direction reverses at every cutoff: later specimens sat in more uncovered climate than the existing reference collection.


![Fig4](figs/Fig4_holdout.png)

**Fig. 4.** Retrospective comparison by collection year. (a) Mean uncovered fraction of later-collected specimens against both nulls; the shaded band marks cutoffs whose training library held fewer than 150 specimens. (b) The same result as a ratio to the recorded-range null, with training library size printed below each point. Values below one indicate later-collected records falling in better-covered climate than a resample of the recorded range. Cutoffs share records and are not independent.

Collecting that supplied later reference specimens has therefore concentrated in climates already better covered than the recorded range at large, while still extending beyond the climates where reference specimens had previously been drawn. The cutoffs share records and are not independent, and records from the same study or locality were not resampled as blocks.

### 3.5. *Oligonychus*

The *Oligonychus* analysis used 2,976 records: 2,798 occurrence records and 178 specimens resolved from 321 rows, 44.6% of which were redundant. Of the 178 specimens, 174 carried a matched verdict, 97.8%; of those, 141 agreed with an anchor and 33 recorded that no anchor was available, so 79.2% of specimens have positive anchor support rather than 97.8%. None fell in a shared cluster. Mean verifiability was 0.844, giving 150.3 effective specimens.

The deficit index was **0.703** counted and **0.745** weighted, a relative widening of 5.9%, with bounds of 0.737 and 0.748. A coverage score below 0.10 covered 1,704 of 2,798 occurrence records, 60.9%, against 9.0% in *Tetranychus*. No cutoff reached the 150-specimen adequacy threshold, so no retrospective comparison is reported for this genus.

The two genera differ in an instructive way. *Oligonychus* is the library with more complete verdict matching and a higher mean evidential weight — 97.8% of its specimens carry a matched verdict and mean verifiability is 0.844, against 94.8% and 0.473 in *Tetranychus* — and the weighting accordingly changes its index less. It is nonetheless far more poorly covered, because it is a smaller library spread across climate space in which three-fifths of occurrence records score below 0.10.

## 4. Discussion

### 4.1. Records, specimens and the cost of not distinguishing them

Resolution to specimens was not a technicality. In *Tetranychus*, 25.3% of rows were additional representations of canonical reference records already present, and the verdict-match rate rose from 70.6% of rows to 94.8% of specimens once identifiers were reconciled. In *Oligonychus* the redundancy was 44.6%. Without resolution, the library would appear larger than it is and much of its existing taxonomic evidence would appear to be missing.

This bears beyond the present analysis. BOLD deposits to GenBank, so any assessment that pools the two databases and keys on source-specific identifiers will inflate its denominator and may mistake an identifier mismatch for absent evidence. Where a crosswalk exists it should be used; where none exists, resulting counts should be reported as records rather than specimens. The same caution applies to our own units: we call them specimens, but establishing that each is one physical mite would require voucher provenance the deposited metadata does not carry.

### 4.2. What the weighting changes, and what it does not

Weighting by verifiability widened the measured gap in *Tetranychus* COI coverage by 28.6%, and the weighted index exceeded the counted index at every coverage scale and under every treatment of records without matched verdicts tested. It is not a finding that public labels are wrong. In this dataset 90.9% of checkable labels were not contradicted, broadly consistent with the general picture for metazoan mitochondrial sequences (Leray et al., 2019) and with the view that reported identification failures often originate in sequence handling rather than in the libraries (Pentinsaari et al., 2020).

The uniform non-contradiction scenario makes the boundary explicit. When every weight is replaced by the bare rate at which labels escape contradiction, the index falls to 0.249 against a counted 0.237 — most of the effect disappears. A framework claiming to measure accuracy would collapse at that point. What survives is the narrower claim: a library user cannot distinguish an unverifiable label that is correct from one that is not, so unverifiable records cannot carry the weight of a reference.

Two cautions attach to that scenario. The rate is estimated only from records whose cluster contains an anchored specimen, and those clusters are the better-studied part of the library, so carrying the rate over to unanchored records is optimistic rather than neutral. And it is a non-contradiction rate at cluster level, not a species-level accuracy: 829 of the 1,563 uncontradicted labels sit in clusters that do not resolve to one species. It should not be compared directly with a genus-level error estimate, nor should the figure of close to 30% unreliable accessions reported by de Mendonça et al. (2011) for a different and much smaller set be transferred to the present library.

### 4.3. Verifiability as a coverage criterion

Existing gap analyses count species with at least one, or at least five, barcode records (Weigand et al., 2019). The threshold convention is sound and we have adopted it, but it is applied to records whose labels are taken as given. For taxa where identification is contested this produces a coverage figure that inherits an unrecorded uncertainty. In *Tetranychus* the barcode region does not separate *T. urticae* from several congeners (Li et al., 2010) and one nominal species has since been synonymised with another (Auger et al., 2013).

Making verifiability an explicit weight carries the taxonomic effort forward into the coverage statistic rather than leaving it in a separate quality assessment. Assessments that quantify library ambiguity and assessments that quantify library coverage have run on parallel tracks — Teixeira et al. (2025), for instance, report both a coverage gap and a 43% rate of BIN-level ambiguity, as separate findings. The present framework combines them into one quantity.

### 4.4. What the index can and cannot support

The environmental-space geometry is not new. Broennimann et al. (2012) established kernel smoothing of occurrence densities in gridded environmental space, with explicit attention to sampling-effort independence. Our contribution is the resolution of records to specimens, the evidential weighting inside that geometry, and the closed-form ranking of candidate localities that follows.

Three limits on interpretation follow from how coverage is computed. Coverage is computed within a genus, so a *T. ludeni* reference contributes where *T. urticae* occurs; the index describes the climatic representation of a genus-wide verified collection and cannot establish species-specific reference support. Extending it to species, or to explicitly defined unresolved complexes, requires harmonised accepted names and is the natural next development. PC1 and PC2 account for 72.7% of bioclimatic variance in *Tetranychus*, so 27.3% lies outside the plane evaluated. And the priority ranking compares alternative single additions against one fixed starting library; it does not constitute an optimal plan for collecting all twenty cells.

### 4.5. The retrospective comparison

Later-collected records sat in significantly less uncovered climate than a resample of the recorded range at every cutoff with an adequate training library. The pattern is that collecting which supplied reference specimens has concentrated where the library was already comparatively strong rather than in the climates it lacked. We have not tested why. Proximity to institutions, concentration in well-studied agricultural regions and existing collaborations are plausible explanations, and each would need its own test.

The comparison is retrospective and narrow in scope. It splits on collection date, so it describes where specimens were gathered, not when they were sequenced or deposited, and it cannot support a claim about how sequencing effort has changed. Cutoffs share records, so the five results are dependent observations of one pattern. The representation and anchor information are fitted on the full dated set, so this is not a prospective validation. And the comparison against the already-sequenced null runs the other way: later specimens did extend beyond the climates from which reference specimens had previously been drawn, just not as far as the recorded range itself extends.

### 4.6. Limitations

WorldClim 2.1 is a 1970–2000 climatology (Fick and Hijmans, 2017), and the records span decades; a specimen from 1960 and one from 2020 receive the same values. At 10 arc-minute resolution a protected-cultivation site and open field 15 km away are indistinguishable, which matters for pests much of whose distribution lies under glass and plastic. The 50 km coordinate-uncertainty limit is coarse relative to that pixel and its effect has not been tested.

Marker eligibility rests on deposited gene or marker fields rather than on alignment to a defined barcode region, so records naming COI without a region were accepted; a confirmed-region-only analysis would establish whether this matters. Coordinate disagreements between two representations of one specimen were resolved by a stated preference rule; four *Tetranychus* groups and 39 *Oligonychus* groups change their climate vector depending on the choice, and the effect of that choice on the index and the priority ranking has not been quantified. We have not tested whether gaps track accessibility and research capacity rather than climate; modelling the deficit against travel time and per-country research output alongside climate would separate those explanations and is the most important extension of this work.

## 5. Conclusions

Public COI reference libraries are routinely assessed by counting records, which presumes both that each record is a distinct specimen and that its species label is correct. Neither presumption held here. A quarter of *Tetranychus* rows and nearly half of *Oligonychus* rows were a second identifier for material already present, and resolving them raised the *Tetranychus* verdict-match rate from 70.6% to 94.8%. Weighting the resulting specimens by verifiability then raised the climate-space coverage deficit from 0.237 to 0.305, a widening that persisted at every coverage scale tested and under both extreme treatments of the 38 reference records without matched verdicts. The smaller *Oligonychus* library had a higher mean evidential weight, but its deficit index rose from 0.703 to 0.745, and 60.9% of its occurrence records sat in climate space with a coverage score below 0.10. A retrospective split on collection year placed later-collected records in significantly better-covered climate than the recorded range at large. The quantity measured is verifiable coverage within a genus, not species-level accuracy and not species-specific reference support: 90.9% of checkable labels were not contradicted, but most of those sit in clusters that do not resolve to one species, and a library user cannot tell which.

---

## Declarations

**Funding.** No funds, grants or other support were received for conducting this study.

**Competing interests.** The authors declare no competing interests.

**CRediT author statement.** **V. K. Chaudhari:** Conceptualization, Methodology, Software, Formal analysis, Data curation, Writing – original draft, Writing – review & editing, Visualization. **Abhishek Shukla:** Conceptualization, Validation, Writing – review & editing, Supervision.

**Data availability.** The analysis code, derived data tables and all results are available at github.com/vkcvaibhav-eng/reference-coverage and will be archived with a citable DOI at acceptance; the release citation will be inserted at that stage. Occurrence data are the GBIF downloads doi:10.15468/dl.xxfcra, doi:10.15468/dl.aejp4w and doi:10.15468/dl.5mckab. Identification verdicts, crosswalks and sequence data are from doi:10.5281/zenodo.23808136.

## Declaration of generative AI in the writing process

The analysis code, figures and a first draft of this manuscript were produced with the assistance of a large language model (Claude, Anthropic). We designed the study, selected and verified the sources, checked the analysis and its outputs, and revised the text. We are responsible for the content.

## References

Auger, P., Migeon, A., Ueckermann, E.A., Tiedt, L., Navajas, M., 2013. Evidence for synonymy between *Tetranychus urticae* and *Tetranychus cinnabarinus* (Acari, Prostigmata, Tetranychidae): review and new data. Acarologia 53, 383–415.

Broennimann, O., Fitzpatrick, M.C., Pearman, P.B., Petitpierre, B., Pellissier, L., Yoccoz, N.G., Thuiller, W., Fortin, M.-J., Randin, C., Zimmermann, N.E., Graham, C.H., Guisan, A., 2012. Measuring ecological niche overlap from occurrence and spatial environmental data. Global Ecology and Biogeography 21, 481–497. https://doi.org/10.1111/j.1466-8238.2011.00698.x

Carné, A., Vieites, D.R., van den Burg, M.P., 2025. In vouchers we (hope to) trust: unveiling hidden errors in GenBank's tetrapod taxonomic foundations. Molecular Ecology 34, e17812. https://doi.org/10.1111/mec.17812

de Mendonça, R.S., Navia, D., Diniz, I.R., Auger, P., Navajas, M., 2011. A critical review on some closely related species of *Tetranychus* sensu stricto (Acari: Tetranychidae) in the public DNA sequences databases. Experimental and Applied Acarology 55, 1–23. https://doi.org/10.1007/s10493-011-9453-5

Fick, S.E., Hijmans, R.J., 2017. WorldClim 2: new 1-km spatial resolution climate surfaces for global land areas. International Journal of Climatology 37, 4302–4315. https://doi.org/10.1002/joc.5086

Leray, M., Knowlton, N., Ho, S.-L., Nguyen, B.N., Machida, R.J., 2019. GenBank is a reliable resource for 21st century biodiversity research. Proceedings of the National Academy of Sciences 116, 22651–22656. https://doi.org/10.1073/pnas.1911714116

Li, G.-Q., Xue, X.-F., Zhang, K.-J., Hong, X.-Y., 2010. Identification and molecular phylogeny of agriculturally important spider mites (Acari: Tetranychidae) based on mitochondrial and nuclear ribosomal DNA sequences, with an emphasis on *Tetranychus*. Zootaxa 2647, 1–15.

Pentinsaari, M., Ratnasingham, S., Miller, S.E., Hebert, P.D.N., 2020. BOLD and GenBank revisited — do identification errors arise in the lab or in the sequence libraries? PLoS ONE 15, e0231814. https://doi.org/10.1371/journal.pone.0231814

Ros, V.I.D., Breeuwer, J.A.J., 2007. Spider mite (Acari: Tetranychidae) mitochondrial COI phylogeny reviewed: host plant relationships, phylogeography, reproductive parasites and barcoding. Experimental and Applied Acarology 42, 239–262. https://doi.org/10.1007/s10493-007-9092-z

Teixeira, M.A.L., et al., 2025. Gaps and data ambiguities in DNA reference libraries: a limiting factor for molecular-based biodiversity assessments using annelids as a case study. Ecology and Evolution 15, e71544. https://doi.org/10.1002/ece3.71544

Weigand, H., Beermann, A.J., Čiampor, F., Costa, F.O., Csabai, Z., Duarte, S., Geiger, M.F., Grabowski, M., Rimet, F., Rulik, B., Strand, M., Szucsich, N., Weigand, A.M., Willassen, E., Wyler, S.A., Bouchez, A., Borja, A., Čiamporová-Zaťovičová, Z., Ferreira, S., Dijkstra, K.-D.B., Eisendle, U., Freyhof, J., Gadawski, P., Graf, W., Haegerbaeumer, A., van der Hoorn, B.B., Japoshvili, B., Keresztes, L., Keskin, E., Leese, F., Macher, J.N., Mamos, T., Paz, G., Pešić, V., Pfannkuchen, D.M., Pfannkuchen, M.A., Price, B.W., Rinkevich, B., Teixeira, M.A.L., Várbíró, G., Ekrem, T., 2019. DNA barcode reference libraries for the monitoring of aquatic biota in Europe: gap-analysis and recommendations for future work. Science of the Total Environment 678, 499–524. https://doi.org/10.1016/j.scitotenv.2019.04.247

Zizka, A., Silvestro, D., Andermann, T., Azevedo, J., Duarte Ritter, C., Edler, D., Farooq, H., Herdean, A., Ariza, M., Scharn, R., Svanteson, S., Wengström, N., Zizka, V., Antonelli, A., 2019. CoordinateCleaner: standardized cleaning of occurrence records from biological collection databases. Methods in Ecology and Evolution 10, 744–751. https://doi.org/10.1111/2041-210X.13152
