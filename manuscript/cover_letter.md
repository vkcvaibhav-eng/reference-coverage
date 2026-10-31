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

Published gap analyses of DNA barcode reference libraries measure coverage by counting records that bear a species name. The count treats that name as given. For taxa where morphological identification is difficult, the assumption is known to fail — in *Tetranychus* sensu stricto, close to 30% of the public sequences examined in one critical review were judged misidentified or dubious, and the barcode region does not separate several of the species involved. A coverage figure built by counting such records inherits that error rate without recording it.

The manuscript sets out a framework that measures instead how much of a taxon's recorded climate space is backed by reference sequences whose species label can be verified against a morphologically identified specimen. Each sequence record carries a verifiability weight; for a record whose sequence cluster holds k anchored species the weight is derived as 1/k rather than chosen. Coverage is then a saturating function of the effective number of verifiable sequences in a neighbourhood of climate space, and the deficit index is the share of recorded occurrence lying outside it. The workflow also returns a ranked list of localities where one further reference specimen would most reduce the shortfall, and a temporal hold-out that tests the surface against sequencing that happened subsequently.

Applied to *Tetranychus*, 974 georeferenced COI records resolved to 405.6 effective sequences and the deficit index rose from 0.218 to 0.282, a relative widening of 29.2% that held between 28.1% and 30.8% across coverage thresholds from k = 1 to k = 20.

Three points about how we have framed the result may interest the reviewers.

First, the paper does not claim that public labels are wrong. In our data 91.2% of checkable labels stood. We state this in the Abstract and report, as an explicit floor, the index obtained when every weight is replaced by the bare probability that a label is correct: under that substitution the effect nearly vanishes. What survives is the narrower and we think defensible claim that a library user cannot distinguish a correct unverifiable label from an incorrect one, so unverifiable records cannot carry the weight of a reference.

Second, we designed the temporal hold-out expecting it to show that subsequent sequencing had been drawn towards the gaps, and it did not. At every cutoff with an adequate training library, post-cutoff sequence records fell in less uncovered climate than draws from the recorded range. We report this as it came out and set out both readings without choosing between them.

Third, the framework requires an independent per-record verification as input. We demonstrate this on a second genus, *Oligonychus*, for which no such verification exists, and report it unweighted with that limitation stated, rather than presenting a weighted figure that would carry no verification information.

All code, derived data tables and results are archived with a citable DOI, the three GBIF occurrence downloads carry their own DOIs, and every script takes a taxon name as an argument, so the analysis reproduces for any group that has a comparable verification.

The manuscript is original, has not been published elsewhere, and is not under consideration by another journal. Both authors have approved the submission. We have no competing interests to declare, and no funding was received for this work. Our use of a generative AI tool in preparing the code and a first draft is declared in the manuscript.

Thank you for considering this work.

Yours sincerely,

V. K. Chaudhari
On behalf of both authors
