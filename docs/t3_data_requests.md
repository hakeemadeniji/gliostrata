# T3 data requests: three emails to send

**To send by 2026-12-01.** Chase non-replies once on 2026-12-15. The registered decision date
is **2027-01-08**: if no qualifying source is secured by then, that T3 row is declared NOT
EVALUABLE (pre-registration item 8), which is a worse outcome than a FAIL because it reports
nothing at all.

Replace `[YOUR NAME]` and `[YOUR AFFILIATION]`. Your email address is in the From field, so
none of these drafts contains it.

**Checked 2026-11-28, and the answer changes email A.** Table S3B is *not* the data this
needs. The paper's own text cites it for "the median time to exceed 2 cm3 ranging from 20 to
189 days" — **one summary value per line**, where the registered criterion needs volumes at
three or more time points per line. Whether the supplement files are downloadable could not
be confirmed either: PMC lists eight files for this manuscript but serves none of them to an
automated request, and the journal's own page returns 403. They are very likely downloadable
in a browser, and worth opening — but on the published text, what is in them is a
time-to-endpoint summary, not serial measurements.

So email A asks for the **underlying serial measurements**, which are probably unpublished,
and says plainly why the published summary will not do. That is a larger ask than pointing at
a table, and it is the honest one.

## What makes a source usable, so you can answer a follow-up question

The registered rule (item 3(b), as amended 2026-11-27) needs, for the tumour-growth row:

- longitudinal volumes — at least three time points per model;
- a control (untreated or vehicle) arm;
- a subtype call or RNA-seq for every model;
- **at least ten models, spanning at least three of the four Verhaak subtypes, with at least
  two models per represented subtype.**

The three-subtype wording matters: modern panels classify with the Wang 2017 scheme, which
has no Neural class, so a three-class panel is expected and is explicitly acceptable. Say so
if asked — it is the difference between a custodian thinking they cannot help and knowing they
can.

Orthotopic implantation is **not** required for this row, so subcutaneous panels qualify.

---

## Email A — Mayo Clinic Brain Tumor PDX National Resource

> **Subject:** Control-arm growth curves from the GBM PDX panel — data request
>
> Dear Dr Sarkaria and colleagues,
>
> I am [YOUR NAME], [YOUR AFFILIATION]. I am running a pre-registered study that tests
> whether a tumour-growth model can reproduce measured PDX growth, and the Mayo panel is the
> best-matched resource I have found.
>
> Vaubel et al. (*Clin Cancer Res* 2020) states that models and genomic data are openly
> available. I have read the paper and its supplementary tables, and I think what I need is
> not among them: Table S3B reports the median time to exceed 2 cm3 per line, whereas the
> analysis needs the serial measurements underlying those curves. Specifically:
>
> - untreated or vehicle control tumour volumes **at each measurement time**, per animal or as
>   a mean per time point, with at least three time points per line — the raw caliper or
>   imaging series rather than a derived time-to-endpoint;
> - the per-line subtype call, or the matching RNA-seq identifier (PRJNA548556) so I can
>   derive it myself.
>
> Ten lines would be sufficient, provided they cover at least three subtypes with two or more
> lines each. A three-class panel under the Wang 2017 scheme is fine — I am not expecting a
> Neural group. If the serial volumes exist only in a lab notebook or an internal database
> rather than in an analysis-ready file, that is still worth knowing: I would rather adapt to
> the format you have than ask you to build one.
>
> The data would be used only to compute a normalised error of a growth model against your
> curves. It would not be redistributed, and the resource would be cited. The analysis plan,
> including the exact criterion above, was registered before any data was sought, and I am
> glad to share it.
>
> Is this something the resource can share, and if so what agreement is needed? If a reply by
> around 19 December is possible it would help, as the study's registered decision point is
> 8 January, after which this analysis is recorded as not evaluable.
>
> With thanks,
> [YOUR NAME]
> [YOUR AFFILIATION]

---

## Email B — Alcañiz et al., Charité / EPO Berlin

> **Subject:** Vehicle-control tumour volumes from your GBM PDX panel (*Front Oncol* 2023)
>
> Dear Dr Alcañiz and colleagues,
>
> I am [YOUR NAME], [YOUR AFFILIATION]. I am running a pre-registered study testing a
> tumour-growth model against measured PDX growth, and your panel of 26 GBM PDX models
> (*Front Oncol* 2023;13:1129627) fits the criteria closely.
>
> The supplement holds the mutation table and the figures, but not the underlying volumes.
> Could you share:
>
> - vehicle-control tumour volumes per measurement time, per animal or as a mean, for as many
>   of the models as is convenient;
> - the per-model subtype assignment already reported in the paper.
>
> Ten models covering three or more subtypes, with two or more per subtype, would be enough.
> Subcutaneous implantation is not a problem for this analysis.
>
> The data would be used only for a calibration check, would not be redistributed, and your
> work would be credited. The analysis plan was registered before the data was sought.
>
> If a reply by around 19 December is possible it would help, as the registered decision point
> is 8 January. If the volumes are not available in a shareable form, a one-line reply saying
> so is genuinely useful — it lets me record the reason rather than leave it open.
>
> With thanks,
> [YOUR NAME]
> [YOUR AFFILIATION]

---

## Email C — in vivo stem-cell migration data

Send to the corresponding author of Rust et al. 2022, copying the Aboody laboratory.

> **Subject:** Quantitative in vivo NSC migration data in glioma — small data request
>
> Dear Dr [corresponding author],
>
> I am [YOUR NAME], [YOUR AFFILIATION]. I am running a pre-registered study of therapeutic
> stem-cell migration toward glioma, and I need measured in vivo data to test a simulation
> against. Rust et al. (2022) notes that raw data are available on request.
>
> What I need is modest:
>
> - per-animal migration distance, or the fraction of delivered cells reaching the tumour,
> - at two or more time points,
> - from MRI, bioluminescence or histology — any of these works.
>
> Aggregate figures published as a single mean are unfortunately not usable, because the
> analysis correlates per-animal values. Even a small dataset, from a handful of animals,
> would be enough to make the comparison possible rather than impossible.
>
> The data would be used only for this correlation, would not be redistributed, and would be
> acknowledged. If you do not have data in this form but know of a group that does, a pointer
> would be just as helpful.
>
> If a reply by around 19 December is possible it would help; the registered decision point is
> 8 January, after which this analysis is recorded as not evaluable.
>
> With thanks,
> [YOUR NAME]
> [YOUR AFFILIATION]

---

## After sending

Record the date each was sent in `docs/t3_data_search.md` §4, and any reply, including a
refusal. Under item 8 both T3 rows are decided independently on 2027-01-08, so a reply on one
row does not hold up the other.
