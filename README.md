# Code accompanying the thesis "Shape Optimisation and Robust Solvers for Incompressible Flow"


This repository contains the code for the shape optimisation examples shown in my thesis.

To obtain Firedrake and the other dependencies, follow the installation instructions below


    # Install firedrake with necessary extra packages
    curl -O https://raw.githubusercontent.com/firedrakeproject/firedrake/master/scripts/firedrake-install
    python3 firedrake-install --opencascade --install pyadjoint --doi 10.5281/zenodo.3369183
    source firedrake/bin/activate
    pip install --no-cache-dir roltrilinos rol

    # Clone and install the Navier-Stokes solver
    git clone https://github.com/florianwechsung/alfi.git
    pip install -e alfi/
    
    # Clone the shape optimisation code and install Fireshape
    git clone --recurse-submodules https://github.com/florianwechsung/ThesisNumerics.git
    cd ThesisNumerics/
    pip install -e fireshape/

