# RealPhy documentation

`Realphy_manual.md` is the source of the RealPhy manual. `Realphy_manual.pdf` is generated from it.

## Building the PDF

Requirements: [pandoc](https://pandoc.org/) and XeLaTeX (TeX Live, MacTeX or BasicTeX) with the packages `fancyhdr`, `titlesec`, `setspace`, `etoolbox`, `xurl`, `xcolor`, `float`, `framed`, `fvextra`, `selnolig`, `unicode-math` and `lm`. The fonts Times New Roman and Courier New are used; the title uses Palatino if available.

On macOS:

```
brew install pandoc
brew install --cask basictex
sudo tlmgr install titlesec xurl framed fvextra selnolig lm lm-math unicode-math
```

Then run

```
./build_pdf.sh
```

in this folder. The result is `Realphy_manual.pdf`.

## Files

| File | Purpose |
|---|---|
| `Realphy_manual.md` | Manual text (Markdown) |
| `img/` | Figures used by the manual |
| `header.tex` | LaTeX settings (spacing, headers, title) used when building the PDF |
| `build_pdf.sh` | Builds the PDF |
