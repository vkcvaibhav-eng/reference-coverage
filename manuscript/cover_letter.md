V. K. Chaudhari
Department of Entomology
N.M. College of Agriculture
Navsari Agricultural University
Navsari 396 450, Gujarat, India
vkcvaibhav@gmail.com

The Editor
Ecological Informatics

Dear Editor,

We submit for your consideration the manuscript "Verifiable coverage of public COI reference libraries in climate space: a framework applied to spider mites (Acari: Tetranychidae)".

Published gap analyses of DNA barcode reference libraries measure coverage by counting records bearing a species name. The count treats two things as given: that each record is a distinct specimen, and that its label is correct. Neither held in our data. GenBank and BOLD hold the same mite under different identifiers, so 974 georeferenced Tetranychus rows resolved to 728 specimens; and the barcode region does not separate several congeners within Tetranychus sensu stricto.

The manuscript sets out a framework that measures instead how much of a taxon's recorded climate space is backed by reference sequences whose species label can be verified against a morphologically identified specimen. Each sequence record carries a verifiability weight; for a record whose sequence cluster holds k anchored species the weight is derived as 1/k rather than chosen. Coverage is then a saturating function of the effective number of verifiable sequences in a neighbourhood of climate space, and the deficit index is the share of recorded occurrence lying outside it. The workflow also returns a ranked list of localities where one further reference specimen would most reduce the shortfall, and a temporal hold-out that tests the surface against sequencing that happened subsequently.

Applied to Tetranychus, 974 rows resolved to 728 specimens, a quarter of the rows being redundant, and resolution raised the verdict-match rate from 70.6% to 94.8%. Weighting the survivors left 344.6 effective specimens and raised the deficit index from 0.237 to 0.305, a relative widening of 28.6%. Oligonychus, 44.6% of whose rows were redundant and 97.8% of whose specimens are verified, returned 0.703 and 0.745 against a much thinner library.

Three points about how we have framed the result may interest the reviewers.

First, the paper does not claim that public labels are wrong. In our data 90.9% of checkable labels were not contradicted, and we report as a sensitivity scenario the index obtained when every weight is replaced by that rate: under that substitution most of the effect disappears. What survives is the narrower claim that a library user cannot distinguish an unverifiable label that is correct from one that is not. We are careful throughout to call this a non-contradiction rate rather than an accuracy, since most uncontradicted labels sit in clusters that do not resolve to a single species.

Second, our retrospective comparison splits on collection year, not on deposition date, and we are explicit that it describes where specimens were gathered rather than how sequencing effort changed. Later-collected specimens fell in significantly better-covered climate than a resample of the recorded range at every adequate cutoff. We designed the test expecting the opposite and report it as it came out.

Third, coverage is computed within a genus, so a congeneric reference contributes to the score where another species occurs. We state plainly that the index describes the climatic representation of a genus-wide verified collection and does not establish species-specific reference support, and we draw no species-level conclusions from it.

This analysis was audited independently twice during preparation. Record-joining, smoothing and archive-loading defects found in those audits are corrected, and every numerical statement in the manuscript is now generated from the released result tables rather than transcribed. All code, derived data and results are archived with a citable DOI, the three GBIF downloads carry their own DOIs, and every script takes a taxon name as an argument.

The manuscript is original, has not been published elsewhere, and is not under consideration by another journal. Both authors have approved the submission. We have no competing interests to declare, and no funding was received for this work. Our use of a generative AI tool in preparing the code and a first draft is declared in the manuscript.

Thank you for considering this work.

Yours sincerely,

V. K. Chaudhari
On behalf of both authors
