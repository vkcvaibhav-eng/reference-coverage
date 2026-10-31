# Verifiable coverage of public COI reference libraries in climate space: a framework applied to spider mites (Acari: Tetranychidae)

V. K. Chaudhari^a,\*^, Abhishek Shukla^b^

^a^ Department of Entomology, N.M. College of Agriculture, Navsari Agricultural University, Navsari 396 450, Gujarat, India

^b^ Department of Entomology, College of Agriculture, Navsari Agricultural University, Waghai Campus 394 730, Gujarat, India

\* Corresponding author. E-mail: vkcvaibhav@gmail.com. ORCID: 0000-0001-9688-0134

---

## Highlights

- Library coverage is weighted by whether each species label can be verified
- Verifiability raises the Tetranychus coverage deficit from 0.258 to 0.327
- One specimen's exact effect on the deficit index ranks candidate localities
- Specimens held under two identifiers inflate counts: 974 rows are 728 mites
- Later-collected specimens sit in better covered climate than the record at large

---

## Abstract

Gap analyses of DNA barcode reference libraries count records, which presumes that the species label on each record is correct. That presumption fails for taxa where morphological identification is hard and public databases carry known misidentifications. We measure instead the *verifiable* coverage of a reference library across the climate space in which a taxon is recorded. Records are first resolved to specimens, since GenBank and BOLD hold the same mite under different identifiers. Each specimen carries a weight for whether its label can be confirmed against a morphologically identified anchor; a specimen whose cluster holds k anchored species is weighted 1/k, derived rather than chosen. The deficit index is the share of recorded occurrence outside coverage; localities are ranked by the exact fall in that index when one verified specimen is added. Applied to *Tetranychus*, 974 georeferenced COI rows resolved to 728 specimens, of which 690 carried a verdict; weighting left 344.5 effective specimens and raised the deficit index from 0.258 to 0.327, a relative widening of 26.9% that held between 25.9% and 38.8% across scales from k = 1 to k = 20. *Oligonychus*, with 174 of 321 specimens verified, gave 0.789 and 0.834. A retrospective split on collection year placed later-collected specimens in significantly less uncovered climate than the recorded range at every adequate cutoff. The quantity measured is verifiability, not accuracy: 90.9% of checkable labels were not contradicted, but a library user cannot tell which.

**Keywords:** DNA barcoding; COI; reference library; coverage; Tetranychidae; biodiversity informatics

---

## 1. Introduction

Species identification from short standardised DNA fragments depends on a reference library in which sequences carry correct species names. Metabarcoding and environmental DNA have deepened that dependence, because the identity assigned to an unknown sequence is inherited from its closest reference match. Assessments of whether reference libraries are adequate have accordingly become a recognised genre. Weigand et al. (2019) analysed coverage of 28,000 aquatic species across BOLD and GenBank for European biomonitoring and reported that coverage varies strongly among taxonomic groups and geographic regions, with many groups holding fewer than five barcodes per species. Comparable regional analyses have since been published for annelids of the Arabian Peninsula, which also found that 43% of the relevant Barcode Index Numbers carried taxonomic ambiguities (Teixeira et al., 2025).

These assessments share a measurement convention: a species is covered when the library holds some number of records bearing its name. The convention treats the name as given. Whether that is safe depends on the taxon. For metazoan mitochondrial sequences as a whole it appears largely safe — Leray et al. (2019) clustered 4.7 million GenBank sequences and estimated an error rate below 1% at the genus level. In difficult groups, and at lower taxonomic levels, the picture differs. Pentinsaari et al. (2020) traced reported identification failures in BOLD and GenBank largely to errors in sequence acquisition and interpretation rather than to the libraries themselves, while Carné et al. (2025) found errors in GenBank's tetrapod records at the level of the vouchers underpinning them.

Spider mites are a group where the convention is demonstrably unsafe. Morphological separation within *Tetranychus* sensu stricto requires slide-mounted males and a small and declining pool of specialists (Ros and Breeuwer, 2007). *Tetranychus urticae*, *T. cinnabarinus*, *T. turkestani* and *T. truncatus* are not reliably separated by COI (Li et al., 2010), and the species status of *T. cinnabarinus* has since been resolved as a red form of the polymorphic *T. urticae* (Auger et al., 2013). Against that background, de Mendonça et al. (2011) evaluated 138 COI and 105 ITS2 GenBank accessions deposited under *T. urticae* and related taxa and judged close to 30% misidentified or dubious. A coverage figure built by counting records bearing a name inherits whatever error rate attaches to that name without recording it.

A second, independent problem attends counting. Reference sequences are not drawn at random from the places a taxon occurs, and occurrence records are themselves spatially and environmentally biased. A coverage statement that is taxonomic only — so many species with so many barcodes — says nothing about whether the sequences that exist come from the range of environments in which the taxon is found. For pests of global distribution, whose populations differ in thermal biology and whose identification matters at ports and in protected cultivation, the environmental dimension is the one that determines whether the library will work where it is needed.

Quantifying occurrence in environmental rather than geographic space is itself established practice. Broennimann et al. (2012) set out a framework that applies kernel smoothers to occurrence densities in gridded environmental space, deliberately so that the result does not depend on sampling effort or on an arbitrary resolution, and used it to compare niches between taxa. We adopt that geometry. What we add is not the environmental-space representation but the evidential weighting within it: the question is not where a taxon's occurrences sit in climate space, but how much of that space is backed by reference sequences a user could actually check.

The framework described here therefore asks, for a given taxon, how much of its recorded climate space is covered by reference sequences whose label can be verified against a morphologically identified specimen, and where one further specimen would most reduce the shortfall. The quantity is deliberately verifiability rather than accuracy. We take no position on how often public labels are right; we measure how much of the library a user could check. Spider mites serve as the demonstration because an independent per-record verification against morphologically identified anchors already exists for them, but every step takes a taxon name as an argument and nothing is specific to Tetranychidae.

## 2. Materials and methods

### 2.1. Occurrence records

Occurrence records were obtained from GBIF through the download API, which issues a citable DOI for a frozen snapshot, on 31 October 2026: *Tetranychus urticae* (6,646 records; doi:10.15468/dl.xxfcra), *Tetranychus* at genus level (15,793 records; doi:10.15468/dl.aejp4w) and *Oligonychus* at genus level (4,706 records; doi:10.15468/dl.5mckab).

Coordinates were screened following the categories of problem catalogued by Zizka et al. (2019). Records were removed when they carried a GBIF issue flag that makes a coordinate unusable (`ZERO_COORDINATE`, `COORDINATE_OUT_OF_RANGE`, `COORDINATE_INVALID`, `COUNTRY_COORDINATE_MISMATCH`, `PRESUMED_SWAPPED_COORDINATE`, `PRESUMED_NEGATED_LATITUDE`, `PRESUMED_NEGATED_LONGITUDE`), when latitude equalled longitude, when both were integers, when the stated coordinate uncertainty exceeded 50 km, or when a flag indicated derivation from an administrative centre. These are operational filters on the fields GBIF supplies, not determinations about how a record was made; `COORDINATE_ROUNDED`, in particular, indicates decimal rounding rather than an established administrative centroid, and we treat it conservatively.

Records whose basis was `HUMAN_OBSERVATION` or `MACHINE_OBSERVATION` were held back from the main analysis and written to separate files. The label denotes an observation made by a person, not that the record is photograph-based, and an expert field determination may legitimately carry it; we hold these records back because the basis field alone does not establish that a specimen was examined. For the genus-level *Tetranychus* download this screen retained 9,623 of 15,793 records, removing 3,894 without coordinates, 1,227 with centroid-indicating flags, 133 with fatal issues, 36 with integer coordinates and 26 with excessive uncertainty, and holding back 854 observations.

### 2.2. Reference sequence records and resolution to specimens

GenBank and BOLD records for *Tetranychus* and *Oligonychus* were taken from the archived deposit of a preceding study of these libraries (doi:10.5281/zenodo.23808136), so that the library measured here is the one the identification verdicts describe. Sequence records were georeferenced from the GenBank `lat_lon` qualifier and from BOLD coordinates; BOLD coordinates flagged as centroid-derived were excluded at extraction. Collection year was recovered from the GenBank `collection_date` qualifier, the BOLD collection date fields and the GBIF `year` and `eventDate` fields.

A record is not a specimen. BOLD deposits its sequences in GenBank, so one mite commonly appears twice: once under a BOLD process ID and once under a GenBank accession. Keying on source and record identifier counts such a mite twice and fails to find a verdict filed under its other identifier. We therefore built an alias map from the crosswalk published by the preceding study, which carries accession and process ID side by side, and resolved every row to a canonical specimen. For *Tetranychus* this collapsed 246 rows, reducing 974 rows to 728 specimens; 101 of the collapsed groups disagreed on coordinates between their two representations, and these are listed in the released results rather than silently reconciled. Among specimens sharing a canonical key, the retained row is the one carrying a verdict, then by source preference, then by identifier, so the outcome does not depend on input order.

Only the barcode region can supply reference coverage for a COI library. Marker was therefore treated as a filter: a georeferenced record whose sequence was COI-3P, another gene, or absent was moved to the occurrence side, where it still locates the taxon but does not count as reference coverage. For *Tetranychus* this moved 94 COI-3P, one CYTB and one 18S record, which is why the 9,915 occurrence records comprise 9,819 GBIF records and 96 reclassified sequence records. Marker assignment rests on the deposited gene or marker field; for records whose field names COI without specifying the region, position within COI was not independently confirmed by alignment, and we flag this as a limitation.

### 2.3. Identification verdicts and verifiability weights

Each specimen carries a per-record verdict from the preceding study, in which public records were tested against morphologically identified anchor specimens drawn from taxonomic publications. The verdicts and their weights are given in Table 1. Verdicts were matched through the alias map, so a verdict filed under either identifier is found: 690 of 728 *Tetranychus* specimens and 174 of 321 *Oligonychus* specimens carried one.

**Table 1.** Verifiability weights. Only the shared-cluster weight is derived from the data. The remainder state what an unverifiable record is worth as a reference and are conventions, not estimates of how often such labels are correct.

| Verdict | Weight | Meaning | *Tetranychus* specimens |
|---|---|---|---|
| Agrees with anchor | 1.00 | label concordant with an anchored specimen in the same cluster | 174 |
| Agrees, shared cluster | 1/k | cluster holds k anchored species, so the label fits any one of them | 388 |
| No anchor | 0.25 | no anchored specimen in the cluster; label unverifiable | 118 |
| Unnamed, one suggested | 0.25 | record unnamed, one species suggested by clustering | 2 |
| Unnamed, several possible | 0.10 | record unnamed, several species possible | 1 |
| Contradicted | 0.00 | label contradicted by an anchored specimen | 7 |
| No verdict | 0.25 | not covered by the anchor check | 38 |

The shared-cluster weight is derived rather than chosen. A specimen whose cluster contains k anchored species is equally consistent with any of them, so its label is worth 1/k. In the present data every such specimen sat in a cluster holding exactly two anchored species, giving 0.50 throughout; the rule is stated because k is not a constant in general. We present 1/k as an ambiguity convention based on the candidate count, not as a calibrated posterior probability: cluster membership supplies how many species are consistent with the label, not equal prior support for each.

Concordance with an anchor is also weaker evidence than re-examination of the specimen itself. It establishes that the label is consistent with the taxonomic determination of another specimen in the same sequence cluster, under the clustering threshold, alignment region and anchor-eligibility rules of the preceding study. We treat it as such throughout.

### 2.4. Climate space

WorldClim 2.1 bioclimatic variables BIO1–BIO19 at 10 arc-minute resolution (Fick and Hijmans, 2017) were sampled at the coordinates of every georeferenced record. The 19 variables were standardised and reduced by singular value decomposition to the first two principal components, fitted jointly across occurrence and sequence records so that both sets occupy one space; for *Tetranychus* these accounted for 47.2% and 24.7% of variance, leaving 28.1% outside the plane in which coverage is evaluated. A 120 × 120 grid was laid over the occupied region with a 5% margin.

### 2.5. Coverage and the deficit index

Occurrence and sequence records were summarised as Gaussian densities on the grid, the sequence density weighted by verifiability. Smoothing covariance is held fixed across weight schemes. This matters: a kernel density estimator that sets its bandwidth from the weighted effective sample size will spread credit differently as well as award less of it, so a weighted and an unweighted surface would differ in smoothing as well as in evidence. We compute the covariance once from the unweighted sequence points by Scott's rule and reuse it for every scheme, so that comparisons isolate the effect of weighting.

Writing *A* for the area of one kernel bandwidth, 2π√det(**H**), and *N* for the sum of verifiability weights, the expected number of verifiable specimens within one bandwidth of a point *z* is

λ(*z*) = *N* · *f*(*z*) · *A*

and the coverage score saturates as that number approaches a stated scale *k*:

*C*(*z*) = 1 − exp(−λ(*z*)/*k*)

Because *N* is a sum of weights rather than a count, reducing confidence in labels lowers the score in absolute terms rather than redistributing it. Two qualifications attach to the interpretation. λ is a sum of distance-decaying kernel similarities, not a count of specimens inside a hard boundary; and at λ = *k* the score is 0.632, so *k* is a saturation scale rather than a threshold. The deficit index is the occurrence-weighted mean shortfall,

*D* = ∫ *p*(*z*)(1 − *C*(*z*)) d*z*

where *p* is occurrence density normalised to unit mass, so *D* = 0 is full coverage and *D* = 1 is none. *D* is a smooth occurrence-weighted shortfall score, not a measured fraction of populations lacking a barcode. We set *k* = 5 for the primary analysis, following the threshold at which Weigand et al. (2019) report barcode coverage alongside the single-record threshold, and report *D* across *k* = 1 to 20.

Coverage is computed within a genus. A *T. ludeni* reference therefore contributes to the score at climates where *T. urticae* occurs. The index measures the climatic representation of a genus-wide verified sequence collection against a genus-wide occurrence record; it does not establish that any named species has species-specific reference support, and we draw no species-level conclusions from it.

### 2.6. Sampling priorities

Occurrence records were binned into 2° geographic cells by the floor of their coordinates, so a cell labelled −14, 130 spans 14° to 12° S and 130° to 132° E. For each cell we compute the exact fall in *D* produced by adding one verified specimen of unit weight at that cell's climate centroid. Because adding one specimen raises λ by exactly *A* times a single kernel evaluated at the centroid, the new index is available in closed form and no approximation or suppression constant is involved. Cells are ranked by that reduction. Uncovered mass and mean uncovered fraction are reported alongside it for every cell, since they answer different questions: uncovered mass favours cells holding many occurrence records, mean uncovered fraction favours climatically unrepresented cells whatever their record count.

### 2.7. Retrospective comparison by collection year

To test whether the deficit surface describes something that later collecting has tracked, the analysis was trained on specimens collected up to a cutoff year and compared against specimens collected afterwards. The split is on collection date, not on deposition or publication date: a mite collected in 2012 and sequenced in 2024 is placed in the training set. The comparison is therefore retrospective with respect to collecting, and it cannot speak to how sequencing effort changed over time.

Two nulls were used, each resampling with equal probability, 2,000 draws, the uncovered fractions of a stated reference set: the pre-cutoff occurrence records (*recorded range*) and the pre-cutoff sequence specimens (*already sequenced*). Weighting such draws by a density estimated from the same points would sample the square of that density and would not correspond to either benchmark. Both tails are reported with the standard Monte Carlo correction (*r* + 1)/(*m* + 1); a large upper-tail *p* is a failure to detect over-targeting and is not evidence of under-targeting. Cutoffs from 2006 to 2020 were swept; those whose training library held fewer than 150 specimens are flagged and not interpreted. Cutoffs share records and the observations are therefore dependent, so the sweep is read as a description of stability rather than as independent replication.

### 2.8. Sensitivity

Weight schemes spanned a naive scheme in which every specimen counts fully, strict and lenient variants, an anchored-only scheme, a scheme crediting unverified specimens in full, and a *non-contradiction floor* in which every weight is replaced by the observed rate at which checkable labels were not contradicted. The coverage scale *k* was swept from 1 to 20. Each run also reports the index with unverified specimens credited in full and discarded, bounding the effect of the default weight.

### 2.9. Software and data availability

The workflow is implemented in Python 3.12 using NumPy, SciPy, rasterio and Biopython, with 30 automated tests. All figures are drawn from the released result tables by the plotting script in the repository. Code, derived data and results are archived at Zenodo (DOI to be inserted at acceptance) and the repository is public at github.com/vkcvaibhav-eng/reference-coverage. Every script takes a taxon name as an argument.

## 3. Results

### 3.1. The analysed set

After screening and resolution to specimens, the *Tetranychus* analysis used 10,643 georeferenced records carrying a complete BIO1–BIO19 vector: 9,915 on the occurrence side (9,819 GBIF records and 96 reclassified non-barcode sequence records) and 728 reference specimens. Resolution to specimens collapsed 246 duplicate rows through an alias map of 3,607 identifiers.

Of the 728 specimens, 690 carried an identification verdict and 38 did not. Mean verifiability was 0.473, so 728 specimens amounted to **344.5 effective specimens**. Across all records in the anchor tables that could be checked against an anchored specimen, 1,563 of 1,720 labels were not contradicted and 157 were, a non-contradiction rate of 0.909. This is not a measured species-label accuracy: 1,009 of the uncontradicted labels sit in clusters holding more than one anchored species and are therefore consistent with, but not resolved to, the deposited name.

### 3.2. Coverage deficit

Counting specimens gave a deficit index of **0.258**. Weighting by verifiability gave **0.327**, a widening of 0.069 or **26.9% in relative terms** (Fig. 1). Under the weighted surface, 911 of 9,915 occurrence records — 9.2% — sat in climate space with a coverage score below 0.10.

The direction of the widening was stable across the coverage scale, though its size was not: across *k* = 1, 3, 5, 10 and 20 the relative widening fell from 38.8% to 25.9% as the index rose from 0.097 to 0.473 unweighted and 0.134 to 0.596 weighted (Fig. 3b). The weighted index exceeded the counted index at every scale tested.

![Fig1](figs/Fig1_climate_space.png)

**Fig. 1.** Climate space of *Tetranychus* on the first two principal components of BIO1–BIO19. (a) Occurrence density from 9,915 georeferenced records, of which 9,819 are GBIF occurrences and 96 are sequence records reclassified as non-barcode. (b) Verifiable coverage score from 728 reference specimens weighted by verifiability, at a saturation scale of five verifiable specimens per kernel bandwidth. (c) The resulting coverage deficit, the product of occurrence mass and the shortfall in the coverage score. White contours show occurrence density in all three panels.


The 38 specimens without a verdict make little difference. Crediting them in full gave 0.306 and discarding them gave 0.337, both above the counted 0.258. The non-contradiction floor, in which every weight becomes 0.909, gave 0.270 (Fig. 3a). An anchored-only scheme, retaining 132 effective specimens, gave 0.493.

![Fig3](figs/Fig3_sensitivity.png)

**Fig. 3.** Sensitivity of the deficit index. (a) Weight schemes, with effective specimen counts in parentheses; the dashed line is the counted index of 0.258. (b) Coverage scale *k* swept from 1 to 20, with the relative widening between the counted and weighted indices annotated at each value. The direction of the widening is stable; its magnitude falls from 38.8% to 25.9% as *k* rises.


### 3.3. Where the gaps are

The 9,915 occurrence records fell in 655 two-degree cells. Ranked by the exact reduction in *D* from one verified specimen, the leading cells were on the United States Gulf coast (30–32° N, 92–90° W; Δ*D* = 0.0078 on 8 occurrence records), the Azores (38–40° N, 28–26° W; 0.0077 on 22 records), south-eastern Australia (34–32° S, 152–154° E; 0.0077 on 51 records), a second Azores cell (38–40° N, 30–28° W; 0.0076 on 194 records) and western Japan (34–36° N, 130–132° E; 0.0072 on a single record) (Fig. 2).

![Fig2](figs/Fig2_geographic_deficit.png)

**Fig. 2.** Coverage deficit mapped back to geography. Each point is a 2° cell holding *Tetranychus* occurrence records, plotted at its centre; colour gives the mean uncovered fraction of the records it contains and symbol size scales with record count. Open circles mark the ten cells where one verified specimen would most reduce the deficit index, numbered for the first five. Note that this ranking is not a ranking by uncovered fraction.


The ranking is not a ranking by shortfall. The leading cell holds an uncovered mass of 4.48, while the fourth holds 107.2; what distinguishes them is where their climate sits relative to the existing verified collection, not how much uncovered occurrence they contain. Cells whose climate is already well represented return small reductions however many uncovered records they hold.

### 3.4. Retrospective comparison by collection year

At the 2015 cutoff, 354 training specimens gave a deficit index of 0.523, and the 289 specimens collected afterwards fell at a mean uncovered fraction of 0.407. A resample of the pre-cutoff occurrence records expected 0.547 and a resample of the pre-cutoff sequence specimens expected 0.109.

Across the five cutoffs with a training library of at least 150 specimens — 2014, 2015, 2016, 2018 and 2020 — later-collected specimens fell in **significantly less** uncovered climate than a resample of the recorded range, with lower-tail *p* = 0.0005 at every cutoff and ratios to that null of 0.75, 0.75, 0.78, 0.65 and 0.53 (Fig. 4). Against the already-sequenced null the direction reverses, with later specimens in more uncovered climate than existing sequencing at every cutoff.


![Fig4](figs/Fig4_holdout.png)

**Fig. 4.** Retrospective comparison by collection year. (a) Mean uncovered fraction of later-collected specimens against both nulls; the shaded band marks cutoffs whose training library held fewer than 150 specimens. (b) The same result as a ratio to the recorded-range null, with training library size printed below each point. Values below one indicate later-collected specimens falling in better-covered climate than a resample of the recorded range. Cutoffs share records and are not independent.

Collecting that supplied later reference specimens has therefore concentrated in climates already better covered than the recorded range at large, while still extending beyond the climates where sequencing had previously occurred. The cutoffs share records and are not independent.

### 3.5. *Oligonychus*

The *Oligonychus* analysis used 3,199 records (2,878 occurrence, 321 specimens), of which 174 carried a verdict. Mean verifiability was 0.579, giving 186.0 effective specimens. The deficit index was **0.789** counted and **0.834** weighted, a relative widening of 5.8%, with bounds of 0.811 and 0.855. A coverage score below 0.10 covered 2,153 of 2,878 occurrence records, 74.8%, against 9.2% in *Tetranychus*.

The widening is smaller than in *Tetranychus* because the index is already close to its ceiling: with so little coverage to begin with, further discounting the library changes the shortfall less.

## 4. Discussion

### 4.1. What the weighting changes, and what it does not

Weighting by verifiability widened the measured gap in *Tetranychus* COI coverage by 26.9%, and the direction of that widening held across every coverage scale and every treatment of unverified specimens tested. It is not a finding that public labels are wrong. In this dataset 90.9% of checkable labels were not contradicted, broadly consistent with the general picture for metazoan mitochondrial sequences (Leray et al., 2019) and with the view that reported identification failures often originate in sequence handling rather than in the libraries (Pentinsaari et al., 2020).

The non-contradiction floor makes the boundary explicit. When every weight is replaced by the bare rate at which labels escape contradiction, the index falls to 0.270 against a counted 0.258 — most of the effect disappears. A framework that claimed to measure accuracy would collapse at that point. What survives is the narrower claim: a library user cannot distinguish an unverifiable label that is correct from one that is not, so unverifiable records cannot carry the weight of a reference.

Two cautions attach to the floor. The rate is estimated only from records whose cluster contains an anchored specimen, and those clusters are the better-studied part of the library, so carrying the rate over to unanchored records is optimistic rather than neutral. And it is a non-contradiction rate at cluster level, not a species-level accuracy: 1,009 of the 1,563 uncontradicted labels sit in clusters that do not resolve to one species. It should not be compared directly with a genus-level error estimate, nor should the historical figure of close to 30% unreliable accessions reported by de Mendonça et al. (2011) for a different and smaller set of accessions be transferred to the present library.

### 4.2. Records, specimens and the cost of not distinguishing them

Resolution to specimens was not a technicality. In *Tetranychus*, 246 of 974 rows were second representations of a mite already present, and 248 of the 286 rows that appeared to lack a verdict in fact carried one filed under their other identifier. Had we not resolved them, the library would have appeared a quarter larger than it is and three times less verifiable than it is.

This has a bearing beyond the present analysis. BOLD deposits to GenBank, so any assessment that pools the two databases and keys on source-specific identifiers will inflate its denominator and may mistake an identifier mismatch for missing evidence. Where a crosswalk exists it should be used; where none exists, the resulting counts should be read as records rather than specimens.

### 4.3. Verifiability as a coverage criterion

Existing gap analyses count species with at least one, or at least five, barcode records (Weigand et al., 2019). The threshold convention is sound and we have adopted it, but it is applied to records whose labels are taken as given. For taxa where identification is contested this produces a coverage figure that inherits an unrecorded uncertainty. In *Tetranychus* the barcode region does not separate *T. urticae* from several congeners (Li et al., 2010) and one nominal species has since been synonymised with another (Auger et al., 2013).

Making verifiability an explicit weight carries the taxonomic effort forward into the coverage statistic rather than leaving it in a separate quality assessment. Assessments that quantify library ambiguity and assessments that quantify library coverage have run on parallel tracks — Teixeira et al. (2025), for instance, report both a coverage gap and a 43% rate of BIN-level ambiguity, as separate findings. The present framework combines them into one quantity.

### 4.4. What the index can and cannot support

The environmental-space geometry is not new. Broennimann et al. (2012) established kernel smoothing of occurrence densities in gridded environmental space, with explicit attention to sampling-effort independence, and the niche-overlap literature has developed it since. Our contribution is the evidential weighting inside that geometry and the closed-form ranking of candidate localities that follows from it.

Two limits on interpretation follow from how coverage is computed. First, it is computed within a genus: a *T. ludeni* reference contributes to the score where *T. urticae* occurs. The index therefore describes the climatic representation of a genus-wide verified collection, and cannot establish that any named species has species-specific reference support. Extending the framework to species, or to explicitly defined unresolved complexes, requires harmonised accepted names and is the natural next development. Second, PC1 and PC2 account for 71.9% of bioclimatic variance, so 28.1% lies outside the plane in which coverage is evaluated.

### 4.5. The retrospective comparison

Later-collected specimens sat in significantly less uncovered climate than a resample of the recorded range, at every cutoff with an adequate training library. The plain reading is that collecting which supplied reference specimens has concentrated where the library was already comparatively strong — near institutions, in well-studied agricultural regions, along established collaborations — rather than in the climates the library lacked.

The test is retrospective and its scope is narrow. It splits on collection date, so it describes where specimens were gathered, not when they were sequenced or deposited, and it cannot support a claim about how sequencing effort has changed. Cutoffs share records, so the five results are dependent observations of one pattern rather than five independent tests. And the comparison against the already-sequenced null runs the other way: later specimens did extend beyond the climates where sequencing had previously occurred, just not as far as the recorded range itself extends.

What the comparison does not do is validate the deficit surface as a predictor. We designed it expecting the opposite result and report it as it came out.

### 4.6. *Oligonychus*

*Oligonychus* has barely half its specimens verified and a counted deficit three times that of *Tetranychus*, with three-quarters of its occurrence records in climate scoring below 0.10 coverage. The weighting changes the index less here, not because verification matters less but because the index is near its ceiling: when little is covered, discounting what exists has less room to act. The contrast between the two genera is the clearest demonstration in this study that the framework distinguishes a library that is thin from one that is thin and unverifiable.

### 4.7. Limitations

WorldClim 2.1 is a 1970–2000 climatology (Fick and Hijmans, 2017), and the records span decades; a specimen from 1960 and one from 2020 receive the same values. At 10 arc-minute resolution a protected-cultivation site and open field 15 km away are indistinguishable, which matters for pests much of whose distribution lies under glass and plastic. We have not separated protected from open cultivation, and a land-cover layer would allow it. The 50 km coordinate-uncertainty limit is coarse relative to that pixel and its effect has not been tested.

Marker eligibility rests on deposited gene or marker fields rather than on alignment to a defined barcode region, so records naming COI without a region were accepted. The 101 specimens whose two database representations disagree on coordinates were resolved by a stated preference rule rather than by adjudication. We have not tested whether gaps track accessibility and research capacity rather than climate; modelling the deficit against travel time and per-country research output alongside climate would separate those explanations and is the most important extension of this work. Finally, the default weight of 0.25 for unverified specimens is a stated position, not a measurement; with only 38 such specimens in *Tetranychus* its influence is small, but in *Oligonychus*, with 147, it is not.

## 5. Conclusions

Public COI reference libraries are routinely assessed by counting records, which presumes both that each record is a distinct specimen and that its species label is correct. Neither presumption held here. Resolving *Tetranychus* records to specimens reduced 974 rows to 728 mites and recovered verdicts for 248 rows that had appeared unverified; weighting the survivors by verifiability then raised the climate-space coverage deficit from 0.258 to 0.327. The widening persisted at every coverage scale tested and under both extreme treatments of the 38 remaining unverified specimens. *Oligonychus*, with half its specimens verified, returned 0.789 and 0.834 against a far thinner library. A retrospective split on collection year placed later-collected specimens in significantly better-covered climate than the recorded range at large, suggesting that collecting has followed existing strength rather than existing gaps. The quantity measured is verifiable coverage within a genus, not species-level accuracy and not species-specific reference support: 90.9% of checkable labels were not contradicted, but a library user cannot tell which, and that is why the distinction is worth drawing.

---

## Declarations

**Funding.** No funds, grants or other support were received for conducting this study.

**Competing interests.** The authors declare no competing interests.

**CRediT author statement.** **V. K. Chaudhari:** Conceptualization, Methodology, Software, Formal analysis, Data curation, Writing – original draft, Writing – review & editing, Visualization. **Abhishek Shukla:** Conceptualization, Validation, Writing – review & editing, Supervision.

**Data availability.** Code, derived data and results are archived at Zenodo and available at github.com/vkcvaibhav-eng/reference-coverage. Occurrence data are the GBIF downloads doi:10.15468/dl.xxfcra, doi:10.15468/dl.aejp4w and doi:10.15468/dl.5mckab. Identification verdicts and sequence data are from doi:10.5281/zenodo.23808136.

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
