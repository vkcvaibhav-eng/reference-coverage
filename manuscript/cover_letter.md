V. K. Chaudhari
Department of Entomology
N.M. College of Agriculture
Navsari Agricultural University
Navsari 396 450, Gujarat, India
vkcvaibhav@gmail.com

The Editor
Ecological Informatics

Dear Editor,

We submit for your consideration the manuscript "Evidence-weighted climatic coverage of public COI reference libraries: a framework applied to spider mites".

Published gap analyses of DNA barcode reference libraries measure coverage by counting records bearing a species name. The count treats two things as given: that each record is a distinct specimen, and that its label is correct. Neither held in our data. GenBank and BOLD index the same reference material under different identifiers, so 974 georeferenced Tetranychus rows resolved to 728 specimens; and the barcode region does not separate several congeners within Tetranychus sensu stricto.

The manuscript sets out a framework that measures instead how much of a taxon's recorded climate space is backed by reference sequences whose species label can be verified against a morphologically identified specimen. Records are first resolved to canonical reference records, since GenBank and BOLD index the same material under different identifiers. Each then carries a verifiability weight, with an equal-ambiguity convention of 1/k where a sequence cluster holds k anchored species. Coverage is a saturating function of the local kernel-weighted reference support, and the deficit index is the occurrence-density-weighted mean shortfall in that coverage score. The workflow also returns a ranked list of localities where one further verified reference record would most reduce the shortfall.

Applied to Tetranychus, 974 rows resolved to 728 specimens, a quarter of the rows being redundant, and resolution raised the verdict-match rate from 70.6% to 94.8%. Weighting the survivors left 344.6 effective specimens and raised the deficit index from 0.237 to 0.305, a relative widening of 28.6%. For Oligonychus, 321 database rows resolved to 178 canonical reference records. Of these, 174 (97.8%) had a matched verdict: 141 agreed with an anchor and 33 had no anchor available. Weighting increased the deficit index from 0.703 to 0.745.

Three points about how we have framed the result may interest the reviewers.

First, the paper does not claim that public labels are wrong. In our data 90.9% of checkable labels were not contradicted, and we report as a sensitivity scenario the index obtained when every weight is replaced by that rate: under that substitution most of the effect disappears. What survives is the narrower claim that a library user cannot distinguish an unverifiable label that is correct from one that is not. We are careful throughout to call this a non-contradiction rate rather than an accuracy, since most uncontradicted labels sit in clusters that do not resolve to a single species.

Second, a retrospective comparison by collection year assesses whether later-collected reference records occupy better- or less-covered climate than the recorded-range benchmark. They fell in significantly better-covered climate at every adequate cutoff. We designed the test expecting the opposite and report it as it came out; it describes where material was gathered and says nothing about how sequencing effort changed.

Third, coverage is computed within a genus, so a congeneric reference contributes to the score where another species occurs. We state plainly that the index describes the climatic representation of a genus-wide verified collection and does not establish species-specific reference support, and we draw no species-level conclusions from it.

This analysis was audited independently during preparation, and record-joining, smoothing and archive-loading defects found in those audits are corrected. The analysis code, derived data tables and results are archived at doi:10.5281/zenodo.23841628, which provides automated tests, figure-generation scripts and machine-readable numerical summaries; development continues at https://github.com/vkcvaibhav-eng/reference-coverage. The occurrence downloads and the preceding reference-library dataset also have citable DOIs.

The manuscript is original, has not been published elsewhere, and is not under consideration by another journal. Both authors have approved the submission. We have no competing interests to declare, and no funding was received for this work. Our use of a generative AI tool in preparing the code and a first draft is declared in the manuscript.

Thank you for considering this work.

Yours sincerely,

V. K. Chaudhari
On behalf of both authors
