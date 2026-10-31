V. K. Chaudhari
Department of Entomology
N.M. College of Agriculture
Navsari Agricultural University
Navsari 396 450, Gujarat, India
vkcvaibhav@gmail.com

The Editor
Biodiversity and Conservation

Dear Editor,

We submit for your consideration the manuscript "Evidence-weighted climatic coverage of public COI reference libraries: a framework applied to spider mites".

Molecular identification now carries much of biodiversity monitoring, and what a metabarcoding or environmental DNA survey can report is bounded by the reference library behind it. Published gap analyses measure that library by counting records bearing a species name. The count treats two things as given: that each record is a distinct specimen, and that its label is correct. Neither held in our data. GenBank and BOLD index the same reference material under different identifiers, so 974 georeferenced Tetranychus rows resolved to 728 specimens; and the barcode region does not separate several congeners within Tetranychus sensu stricto.

The manuscript sets out a framework that measures instead how much of a taxon's recorded climate space is backed by reference sequences whose species label can be verified against a morphologically identified specimen. Records are first resolved to canonical reference records. Each then carries a verifiability weight, with an equal-ambiguity convention of 1/k where a sequence cluster holds k anchored species. Coverage is a saturating function of the local kernel-weighted reference support, and the deficit index is the occurrence-density-weighted mean shortfall in that coverage score.

What we think will interest this journal's readers is the output that follows. Because the fall in the index from one added specimen has a closed form, candidate localities can be ranked, and the ranking answers a question a count-based gap analysis cannot: not whether more barcodes are needed, but where a collecting effort would buy the most. For Tetranychus the leading cells span 6 to 177 occurrence records, so neither the emptiest cells nor the best-recorded ones win automatically.

Applied to Tetranychus, a quarter of the rows were redundant, and resolution raised the verdict-match rate from 70.6% to 94.8%. Weighting left 344.6 effective specimens and raised the deficit index from 0.237 to 0.305, a relative widening of 28.6%. For Oligonychus, 321 database rows resolved to 178 canonical reference records, of which 174 (97.8%) had a matched verdict: 141 agreed with an anchor and 33 had no anchor available. Weighting raised its index from 0.703 to 0.745, and three-fifths of its occurrence records sit in climate scoring below 0.10 coverage against 9.0% in Tetranychus.

Three points about how we have framed the result may interest the reviewers.

First, the paper does not claim that public labels are wrong. In our data 90.9% of checkable labels were not contradicted, and we report as a sensitivity scenario the index obtained when every weight is replaced by that rate: under that substitution most of the effect disappears. The claim the index does support is narrower, that a library user cannot distinguish an unverifiable label that is correct from one that is not. We call this a non-contradiction rate rather than an accuracy throughout, since most uncontradicted labels sit in clusters that do not resolve to a single species.

Second, a retrospective comparison by collection year assesses whether later-collected reference records occupy better- or less-covered climate than the recorded-range benchmark. They fell in significantly better-covered climate at every adequate cutoff. We designed the test expecting the opposite and report it as it came out; it describes where material was gathered and says nothing about how sequencing effort changed.

Third, coverage is computed within a genus, so a congeneric reference contributes to the score where another species occurs. We state plainly that the index describes the climatic representation of a genus-wide verified collection and does not establish species-specific reference support, and we draw no species-level conclusions from it.

Spider mites are the demonstration taxon rather than the subject. Every script takes a taxon name as an argument, and applying the framework elsewhere requires occurrence and sequence inputs for that taxon and an equivalent per-record verification, not changes to the code.

This analysis was audited independently during preparation, and record-joining, smoothing and archive-loading defects found in those audits are corrected. The analysis code, derived data tables and results are archived at doi:10.5281/zenodo.23841628, which provides automated tests, figure-generation scripts and machine-readable numerical summaries; development continues at https://github.com/vkcvaibhav-eng/reference-coverage. The occurrence downloads and the preceding reference-library dataset also have citable DOIs.

The manuscript is original, has not been published elsewhere, and is not under consideration by another journal. Both authors have approved the submission. We have no competing interests to declare, and no funding was received for this work. Our use of a generative AI tool in preparing the code and a first draft is declared in the manuscript.

Thank you for considering this work.

Yours sincerely,

V. K. Chaudhari
On behalf of both authors
