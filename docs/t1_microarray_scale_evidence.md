# T1 subtype row: value scale of CGGA mRNA_array_301 (documentation evidence)

Compiled 2026-11-27 from public web pages only. No expression file was downloaded or read.
Pre-registration section 2.2 item 9(d) rule: values enter as distributed if the source's
documentation says they are log2 values, otherwise as log2(x + 1), decided from the
documentation before any value is read.

## What the documentation says

Platform and processing (all verbatim):

1. CGGA download page, https://www.cgga.org.cn/download.jsp, dataset mRNA-array_301:
   data type "mRNA microarray"; platform "Agilent Whole Human Genome (Array)"; "Total number
   of samples: 301"; downloads "Expression Data (gene level)" and "Expression Data (probe
   level)". The page gives no normalisation method, no log2 statement and no units.
2. CGGA about page, https://www.cgga.org.cn/about.jsp: "Microarray analysis was performed
   on all samples using the Agilent Whole Human Genome Array"; "Probe intensities were
   normalized using GeneSpring GX 11.0."; "There were 301 samples for mRNA microarray
   detection". It does not mention log2.
3. Zhao Z et al. 2021, Genomics Proteomics Bioinformatics, PMC8498921,
   https://pmc.ncbi.nlm.nih.gov/articles/PMC8498921/ : "cDNA and biotinylated cRNA were
   synthesized and hybridized to the Agilent Whole Human Genome Array"; "the array-generated
   data were analyzed by the Agilent G2565BA Microarray Scanner System"; "GeneSpring GX11.0
   was applied to calculate probe intensity." No normalisation detail, no log2, no GEO
   accession (sequence data go to the Genome Sequence Archive).
4. Yan W et al. 2012, Neuro-Oncology, "Molecular classification of gliomas based on whole
   genome gene expression: a systematic report of 225 samples from the Chinese Glioma
   Cooperative Group", https://pmc.ncbi.nlm.nih.gov/articles/PMC3499016 : "Microarray
   analysis was performed on all 225 samples using the Agilent Whole Human Genome Array";
   "Probe intensities were normalized using GeneSpring GX 11.0." No log2 or other
   transformation is mentioned, and no GEO accession.
5. A third-party meta-analysis, Wang S et al. 2017, Sci Rep, PMC5603565,
   https://pmc.ncbi.nlm.nih.gov/articles/PMC5603565 : "Every dataset was preprocessed using
   log2 transformation". This is the authors' own step applied to every dataset they
   pooled, including CGGA; it is not evidence about how CGGA distributes the values, and it
   points the other way if anything (they log-transformed CGGA themselves).

GEO: no GEO deposit of the 301-sample batch was found. The 225-sample 2012 paper gives no
accession. A search turned up no GEO series for it (GSE16011 is an Affymetrix U133 Plus 2.0
dataset, a different cohort).

## Finding

No CGGA source says the mRNA_array_301 values are log2, and none says they are not. The
documentation names GeneSpring GX 11.0 normalisation only. The pages fetched here were read
through a summarising fetch tool, so these quotes are the tool's extracts of the pages, not
character-checked copies; recheck the wording against the live pages when citing.

## Recommendation

`log2_plus_one`. Confidence: high that this is what the registered rule requires (the
documentation does not state log2, so the "otherwise" branch applies); low to moderate
that the distributed values are in fact linear. GeneSpring workflows often report
log2-scale normalised signal, and Agilent one-colour data are commonly analysed as log2,
so the files may well be log2 already. That is background knowledge I did not verify from a
source, and the rule says it must not decide the scale without documentation.

Consequence if the files are already log2 and log2(x + 1) is applied: the transform is
monotone, so the rank order of genes within a sample is kept, but the z-scored profile is
compressed at the low end, which can shift some Pearson calls. A reviewer may ask for a
sensitivity run on the other scale; that would have to be labelled exploratory and must not
replace the registered scale. Any check of the value range is itself a read of the matrix
and must go through the read guard.

## Open questions

- Whether a CGGA readme inside the downloaded zip states the scale. Reading it is a read of
  the sealed file and is not allowed before the guard call.
- Whether the live download page carries wording not captured by the fetch tool (for
  example a tooltip or a PDF description).
