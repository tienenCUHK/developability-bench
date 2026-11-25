# Building the manuscript

The manuscript uses the official NeurIPS formatting files. `neurips_2026.sty` is included in this directory
(it comes from the official NeurIPS formatting zip; swap in the current year's file when the venue publishes it).

```
latexmk -pdf -bibtex main.tex
```

Tested with TeX Live 2026 (pdfTeX 3.141592653-2.6-1.40.29). Required CTAN packages beyond a minimal scheme:
`booktabs`, `natbib`, `hyperref`, `caption`, `environ`, `lineno`, `psnfss`, `nicefrac`, `units`, `microtype`,
`tools`, `xcolor`, `bm`, `amsmath`, `amssymb`, `graphicx`.

The preamble uses `\usepackage[eandd]{neurips_2026}`, i.e. the anonymous submission mode for the Evaluations &
Datasets track (anonymous author block, line numbers, "Do not distribute" footer). For a non-anonymous local
reading copy, use `\usepackage[preprint,eandd]{neurips_2026}`.

Two checkers are included and do not need a TeX installation:

- `python3 check_latex.py` — environment/brace/math balance, `\ref` targets, graphics paths, `\input` files,
  citation keys against `refs.bib`.
- `python3 verify_against_source.py` — compares every table cell and every number against the markdown source.

The compiled `main.pdf` in this directory is the submission build: main text 9 pages, references and appendices
after, no LaTeX errors, no undefined references.
