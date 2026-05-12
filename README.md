# EpiTables

Interactive summary measures for epidemiology and biostatistics. Runs entirely in the browser — no server, no R, no dependencies beyond a small self-contained stats library.

**Live site:** https://rajsubediresearch.github.io/epitables

---

## Tools
| Page | What it does |
|---|---|
| `index.html` | 2×2 table analysis — OR, RR, PR, RD, NNT, AF, PAF, chi-square, Fisher's exact |
| `diagnostic.html` | Diagnostic accuracy — sensitivity, specificity, PPV, NPV, LR+, LR−, DOR, Youden's index, accuracy |
| `predictive.html` | PPV/NPV prevalence adjustment — sliders for sensitivity, specificity, and population prevalence |
| `mcnemar.html` | McNemar's test for paired proportions — corrected and uncorrected chi-square, difference in positive rates |
| `proportion.html` | Single proportion CI — Wilson score, Wald, and Clopper-Pearson exact methods compared side by side |

---

## Design principles

- **Summary data only** — no individual-level data needed, no file uploads
- **Reactive** — all outputs update instantly as you type or slide
- **Interpretive** — every metric includes a one-liner in plain language using your variable names
- **No backend** — pure HTML + CSS + JS, deployable on GitHub Pages for free
- **Self-contained stats** — `js/stats.js` implements all distributions and CI methods from scratch (no jStat, no R)

---

## Citation

If you use EpiTables in teaching or research, please cite:

Subedi, R. (2026). EpiTables: Interactive summary measures for epidemiology and biostatistics (v1.0.1). Zenodo. https://doi.org/10.5281/zenodo.20135874

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20135874.svg)](https://doi.org/10.5281/zenodo.20135874)

---

## License

MIT
