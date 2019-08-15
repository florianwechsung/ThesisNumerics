# Code accompanying the thesis "Shape Optimisation and Robust Solvers for Incompressible Flow"


This repository contains the code for the shape optimisation examples shown in my thesis.

To obtain Firedrake and the other dependencies, follow the installation instructions below


    curl -O https://raw.githubusercontent.com/firedrakeproject/firedrake/master/scripts/firedrake-install
    python firedrake-install --opencascade --install pyadjoint --doi 10.5281/zenodo.3369183
    pip install --no-cache-dir roltrilinos rol

    git clone --recurse-submodules https://github.com/florianwechsung/ThesisNumerics.git
    cd ThesisNumerics/
    pip install -e fireshape/

