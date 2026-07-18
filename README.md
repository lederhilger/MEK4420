# MEK4420
Notes and assignments based on those of the course [MEK4420 — _Marine Hydrodynamics_](https://www.uio.no/studier/emner/matnat/math/MEK4420/index-eng.html), taught at the University of Oslo.

## On the content
Shared are the notes compiled by Simon Lederhilger, based on the books _Marine Hydrodynamics_ by John Nicholas Newman and _Irrotational Flow of Frictionless Fluids, Mostly of Invariable Density_ by Earle Hesse Kennard.
The content of this repository is relevant to the special curriculum MEKSP100 _The response of wind turbines in waves_, supervised by Tor Anders Nygaard, consisting of four assignments:
* [MA1](https://github.com/lederhilger/MEK4420/blob/main/MA1/latex/main.pdf): _Numerically determining added mass in an unbounded fluid_
* [MA2](https://github.com/lederhilger/MEK4420/blob/main/MA2/latex/main.pdf): _The forces and response of a heaving rectangle_
* [MA3](https://github.com/lederhilger/MEK4420/blob/main/MA3/latex/main.pdf): _The mooring dynamics of floating wind turbines_
* MA4: _Interaction between bodies and irregular sea_

Focus is given to the mathematical theory behind the results, with the numerical results ultimately being secondary.
This theory is founded in complex analysis, here drawing from _Methoden der Komplexen Funktionentheorie_ (originally in Russian: _Методы теории функций комплексного переменного_) by Michael Alekseyevich Lavrentev & Boris Vladimirovich Shabat.
Other sources are found in the [bibliography](https://github.com/lederhilger/MEK4420/blob/main/_latexmk/bibliography.bib).

Further development of the content will follow sporadically.

## Compilation
This project was developed using [TeX Live](https://www.tug.org/texlive/), so an installation of that is recommended.
To compile the project documents (`main.tex`) locally, run the following command from within the directory containing both `main.tex` and `latexmkrc`:

```bash
latexmk main.tex
```

Compiling the figures is done from within the `figures/` directory:

```bash
latexmk -cd path/to/file.tex
```