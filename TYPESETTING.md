# Rebuilding the manuscript

From the repository root:

```sh
latexmk -pdf -interaction=nonstopmode -halt-on-error paper.tex
```

The supplied `paper.pdf` was built with pdfLaTeX on 20 October 2026. The TeX
source is authoritative; typesetting does not run the scientific replay.
