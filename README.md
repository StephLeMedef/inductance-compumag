# Inductor optimization


[![GitHub license](https://img.shields.io/github/license/StephLeMedef/inductance-compumag)](https://github.com/StephLeMedef/inductance-compumag) [![GitHub release](https://img.shields.io/github/release/StephLeMedef/inductance-compumag.svg)](https://github.com/StephLeMedef/inductance-compumag/releases/)

Parametric and free-shape air gap optimization of an inductor using [NGSolve](https://www.ngsolve.org/) and [Mmg](http://www.mmgtools.org/).

## Quickstart

Start by running the Jupyter notebooks in the `notebooks` folder. 

## Contents

- dependencies
    - mmg
      - `mmg2d.dll` : dll file related to remeshing free-shapes with Mmg
- notebooks
    - `FREEFORM.ipynb`, `PARAMETRIC.ipynb`, `REFERENCE.ipynb`, `TAYLOR-TESTS.ipynb` : runable Jupyter notebooks
    - `mmglib.py` : binder to Mmg utilities
    - `shapeOptInductance.py` : module containing useful functions for shape optimization
- `AUTHORS `
- `LICENSE`
- `Pipfile` : description of the environment (items can be installed with `pip`)
- `README.md`

**NB** The code has been tested with success on Windows 11. It should also run on Linux except for the free-form due to the compilation of Mmg library.

## License

Copyright (C) 2025 Stéphane Gaydier (stephane.gaydier@grenoble-inp.org), Itaï Zehavi (itai.zehavi@ens-paris-saclay.fr), Théodore Cherrière (theodore.cherriere@centralesupelec.fr), Peter Gangl (peter.gangl@ricam.oeaw.ac.at).

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.
