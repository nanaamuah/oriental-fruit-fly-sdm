# Climatic suitability for the oriental fruit fly across West Africa, with spatial block validation

A species distribution model for *Bactrocera dorsalis* (Hendel) fitted on public occurrence
records, scored twice: once with conventional random cross-validation, and once with folds
separated in space. The second number is the one that matters.

---

## Read this before reading anything else

This is a model of where climate resembles the climate at places where someone recorded this
insect and uploaded the record. It is not a map of where the insect is, and it is not a pest
risk analysis.

GBIF is a convenience sample. Recording effort follows funded projects, road networks and
institutional interest, so the model's background is partly a map of entomologists. The
target-group background correction used here reduces that problem but does not remove it.

A climate-only model omits host availability, irrigation, trade pathways, natural enemies and
management, all of which decide whether a climatically suitable place is actually infested.
MaxEnt output is a relative suitability index, not a probability of presence, and thresholding
it into a binary suitable/unsuitable map imposes a decision the model does not contain.

The performance drop between random and spatially blocked validation reported below is the
honest number. If you see only the random-fold figure quoted anywhere, it is the wrong one.

Finally, absence of records from a country is not absence of the pest, and that asymmetry falls
hardest on exactly the countries with the least surveillance capacity.

---

## The result

The same MaxEnt model, scored two ways:

| Algorithm | Validation | AUC | TSS | Boyce | AUC drop | TSS drop | Boyce drop |
|---|---|---|---|---|---|---|---|
| MaxEnt | Random | 0.964 ± 0.005 | 0.826 ± 0.014 | 0.984 ± 0.014 | — | — | — |
| MaxEnt | Spatial | 0.790 ± 0.179 | 0.465 ± 0.306 | 0.693 ± 0.437 | 0.174 | 0.361 | 0.291 |
| GLM | Random | 0.947 ± 0.006 | 0.807 ± 0.010 | 0.951 ± 0.024 | — | — | — |
| GLM | Spatial | 0.798 ± 0.167 | 0.507 ± 0.371 | 0.637 ± 0.414 | 0.149 | 0.300 | 0.315 |

Values are the mean and standard deviation across five folds.

Two things to notice. The drop itself, about 0.17 AUC and 0.36 TSS for MaxEnt. And the standard
deviation, which grows by more than an order of magnitude under blocked folds. Under random
folds the five MaxEnt AUCs sit between 0.958 and 0.968. Under blocked folds they run 0.977,
0.553, 0.692, 0.775 and 0.954.

Both algorithms drop by a similar amount, which is the useful part of running two. The
optimism is a property of the validation design, not of MaxEnt.

What that spread means is a separate question from whether it is real, and the report is now
careful about the difference. See "What was resolved" below.

---

## What is in here

```
oriental_fruit_fly/
├── r_scripts/
│   ├── 00_setup.R              environment checks, run once
│   ├── 01_downloads.R          GBIF, WorldClim and EPPO retrieval
│   ├── 02_data_cleaning_qc.R   cleaning, thinning, stopping-rule check
│   ├── 03_analysis.R           predictors, models, both validation schemes
│   └── 04_figures.R            figures, maps and tables
├── data/
│   ├── raw/                    downloaded, never edited (not tracked)
│   └── processed/              script output (tracked, except teph_clean.csv)
├── outputs/
│   ├── figures/                validation plot, response curves
│   ├── maps/                   West Africa and Ghana suitability
│   └── tables/                 everything the report reads
├── qc_plots/                   exploratory cleaning diagnostics (not tracked)
├── report.qmd                  the manuscript, renders to index.html
├── references.bib              bibliography
├── renv.lock                   package versions
├── LICENCE
└── README.md
```

Run the scripts in numerical order. Each writes what the next one reads, so nothing needs to be
held in memory between them.

---

## What is tracked

`.gitignore` is the authority here, and it says the following.

Tracked: `r_scripts/`, `report.qmd`, `references.bib`, `renv.lock`, `LICENCE`, `.Rprofile`,
all of `outputs/`, and all of `data/processed/` apart from one file.

Ignored, with reasons: `data/raw/`, because of source licence conditions and size;
`data/processed/teph_clean.csv`, because it is a large intermediate that
`02_data_cleaning_qc.R` rebuilds in one pass; `qc_plots/`, because those are scoping
diagnostics and not finished outputs; `renv/library/` and `renv/staging/`, restored from
`renv.lock`; and the usual Quarto, R session and OS files.

`data/processed/` being tracked is what allows a reviewer to clone the repository and rerun
`04_figures.R` without refitting the models. The report itself reads only from
`outputs/tables/`, so it renders from a clone without touching `data/processed/` at all.

---

## Data

| Source | What | Licence |
|---|---|---|
| GBIF occurrence download | *B. dorsalis*, 11,757 georeferenced records. DOI [10.15468/dl.jzyuk9](https://doi.org/10.15468/dl.jzyuk9) | Record-level, mixed CC0 / CC BY / CC BY-NC |
| GBIF occurrence download | Tephritidae excluding *B. dorsalis*, 328,343 records, for the target-group background. DOI [10.15468/dl.g9edqf](https://doi.org/10.15468/dl.g9edqf) | As above |
| WorldClim v2.1 | 19 bioclimatic variables at 2.5 arc-minutes | CC BY-SA 4.0 |
| EPPO Global Database | Recorded distribution for DACUDO | EPPO terms, attribute |

Downloading from GBIF through `occ_download()` needs credentials in `.Renviron`
(`GBIF_USER`, `GBIF_PWD`, `GBIF_EMAIL`). That route is used because it mints a citable DOI,
which is what makes the download reproducible by someone else. `occ_search()` does not.

---

## How the presences were reduced from 11,757 to 1,089

| Stage | Records |
|---|---|
| Raw GBIF occurrences | 11,757 |
| After coordinate artefact filtering | 11,709 |
| After removing occurrences without WorldClim coverage | 11,563 |
| After thinning to one occurrence per WorldClim cell | 1,089 |

Three decisions in that table are worth stating openly.

**48 records, not 51.** CoordinateCleaner flagged 35 administrative centroids and 16
biodiversity institutions. Three records were caught by both tests, so the loss is 48.

**Exact-coordinate duplicates were not removed as a separate step.** 856 coordinates carried
repeated records, 10,343 records in total, with one coordinate holding 351. Thinning at
WorldClim cell level handles both exact repeats and distinct coordinates carrying identical
climate values, so a separate duplicate filter would have been redundant.

**Capital-city and sea flags were not used as deletion criteria.** Seven tests were run and
five were used to delete. The capital test flagged 433 records, most of them plausible urban
observations. The sea test proved highly sensitive to coastline resolution. Neither is a
specific enough indicator of error to delete on, so whether a location could enter a climatic
model was decided by whether WorldClim returned a value for it. That is the 146 records in row
three.

The Tephritidae background went through the same five tests, 328,343 down to 324,912, then the
same cell-level thinning. Presences and background must be cleaned identically, or any
difference the model finds between them could be an artefact of the cleaning.

---

## The stopping rule, and what it said

The plan set a checkpoint at the end of day two. If the cleaned African records occupied fewer
than about 25 distinct 100 km cells, the model would be describing a handful of field campaigns
rather than a climate envelope, and the species would be switched to *Spodoptera frugiperda*.

75 cleaned African records fell in 75 distinct 100 km cells. The threshold was cleared and the
project continued with *B. dorsalis*. The check is the last section of `02_data_cleaning_qc.R`
and reruns on demand.

Note what 75 records means in context: 11,757 records globally, 75 of them in Africa after
cleaning, for a pest present in 35 sub-Saharan countries. The model is fitted globally for that
reason and projected onto West Africa, and it is not fitted on African records alone.

---

## Method

Ten of the 19 bioclimatic variables survived VIF screening at a threshold of 10: BIO1, BIO2,
BIO3, BIO8, BIO9, BIO13, BIO14, BIO15, BIO18 and BIO19. Screening was done on the West African
stack rather than globally, because collinearity between bioclimatic variables is regional.

Background points are the recorded locations of other Tephritidae, sampled to 10,000, with any
cell holding a *B. dorsalis* presence excluded. Fruit flies are collected by the same traps and
the same people, so where other Tephritidae were recorded approximates where anyone was
looking. A uniform random background would have handed the model survey effort and let it
report the result as climate.

MaxEnt and a ridge-penalised logistic GLM were both fitted, then both scored under both
schemes. Block size came from the empirical range of spatial autocorrelation in the predictors,
estimated with `cv_spatial_autocor()`, not chosen by eye. `03_analysis.R` writes that number to
`data/processed/block_size_record.csv` and `04_figures.R` carries it into
`outputs/tables/block_size_sensitivity.csv`, so the size actually used is reported and not left
inside the session. The blocked validation was repeated at half and double that size, and the
report tabulates all three.

Within every fold the classification threshold is taken from the training data and applied to
the held-out data. Taking it from the test data would let the threshold adapt to the answer and
would hide the effect being measured.

One caveat about the metrics. The quantity entering TSS as specificity is the proportion of
held-out background points below the threshold, and background points are surveyed localities,
not confirmed absences. These values are comparable between the two validation schemes and are
not comparable with a TSS computed from presence-absence data.

---

## Which predictors, and a caution about them

| Rank | Predictor | Contribution (%) |
|---|---|---|
| 1 | BIO13: Wettest-month precipitation | 50.7 |
| 2 | BIO1: Annual mean temperature | 29.8 |
| 3 | BIO8: Temperature of wettest quarter | 8.9 |
| 4 | BIO3: Isothermality | 4.0 |
| 5 | BIO2: Mean diurnal range | 2.0 |
| 6 | BIO18: Warmest-quarter precipitation | 1.9 |
| 7 | BIO19: Coldest-quarter precipitation | 1.0 |
| 8 | BIO9: Temperature of driest quarter | 0.9 |
| 9 | BIO15: Precipitation seasonality | 0.6 |
| 10 | BIO14: Driest-month precipitation | 0.3 |

Wettest-month precipitation and annual mean temperature together account for four fifths. That
is consistent with what is known about the species, which needs warmth and moist fruiting
seasons. It is also consistent with sampling bias, because people and their fruit trees are
also concentrated in warm, wet places. Percent contribution in MaxEnt is a heuristic that
depends on the order variables entered the fit, so treat the table as a ranking rather than a
partition of explained variance.

---

## Maps and the EPPO overlay

The suitability maps carry the EPPO recorded distribution as country outline colour. EPPO never
entered the model, so it is a genuinely external comparison, but it is a record of reporting
and not a record of presence.

At country level in West Africa it discriminates nothing, because every West African country
carries an EPPO record for this species. The no-record class the overlay was built to reveal is
empty for this region. That is worth stating and not glossing over: the open question in West
Africa is where within a country the risk sits, and a country-level phytosanitary register
cannot answer it.

Saint Helena is in Natural Earth's Western Africa subregion and sits far out in the Atlantic.
It is dropped for plotting only.

---

## What was resolved

Four things were carried as loose ends in earlier versions of this README and of the report.
They are recorded here because the fixes touch tracked outputs.

**The block-size sensitivity test was run and never reported.** Both README and report told the
reader that half and double block sizes had been tested, then pointed at a file under
`data/processed/` without stating the outcome. `03_analysis.R` now also writes the block sizes
themselves, and `04_figures.R` builds `outputs/tables/block_size_sensitivity.csv`, which the
report tabulates. The conclusion holds at all three sizes.

**Fold composition was never reported, and the between-fold spread was over-interpreted.** The
report read the spread in blocked-fold AUC as a map of where the model works and where it
fails. That reading is possible, and so are two others: blocked folds force extrapolation into
climate space the training fold may not cover, and blocked folds hold unequal numbers of
presences, so a small fold returns an unstable AUC for reasons unconnected to the model. Every
CV loop in `03_analysis.R` now records `n_presence` and `n_background` per fold,
`04_figures.R` writes `outputs/tables/fold_composition.csv`, and the report sets out all three
explanations instead of asserting one.

**The report could not render from a clone.** It read `cleaning_log.csv` from
`data/processed/`, while every other table it used came from `outputs/tables/`. All tables the
report reads are now written to `outputs/tables/` by `04_figures.R`, and the report reads
nothing else.

**The renv caveat was false.** An earlier version of this README warned that `renv.lock` did
not list `blockCV`, `modEvA`, `Hmisc`, `lwgeom` or `tidyr`. The committed lockfile pins all
five, at blockCV 4.0-0, modEvA 3.45, Hmisc 5.2-6, lwgeom 0.2-17 and tidyr 1.3.2, and it covers
every package the scripts load. The warning is removed.

---

## Reproducing this

R 4.6.1. Restore the package library with `renv::restore()`, then run `r_scripts/` in order and
render the report:

```r
renv::restore()
source("r_scripts/00_setup.R")
source("r_scripts/01_downloads.R")
source("r_scripts/02_data_cleaning_qc.R")
source("r_scripts/03_analysis.R")
source("r_scripts/04_figures.R")
quarto::quarto_render("report.qmd")
```

`00_setup.R` calls `renv::init()`, which is a one-time step; on a restored project run
`renv::status()` instead. The report is set to `output-file: index.html` with
`embed-resources: true`, so rendering produces one self-contained file and no
`report_files/` directory.

One reproducibility caveat remains. MaxEnt runs on Java through `rJava`, and a version mismatch
does not surface until `MaxEnt()` is first called. `00_setup.R` tests for it up front. If Java
will not cooperate, the ridge GLM is pure R and reproduces the same validation finding on its
own.

Seeds are set before background sampling, fold assignment and every GLM fit. MaxEnt itself is
deterministic given its inputs.

---

## Licence

Code released under the terms in `LICENCE`. Data carry the licences of their sources, listed
above. GBIF downloads must be cited by their DOIs.
