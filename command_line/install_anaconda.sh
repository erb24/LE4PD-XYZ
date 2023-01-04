#!/bin/bash

# install anaconda
ANACONDA="Anaconda3-2022.10-Linux-x86_64.sh"
wget https://repo.continuum.io/archive/${ANACONDA}
bash ${ANACONDA}

# install dependencies
bash<<EOF
pip install notebook
pip install physt
pip install mdtraj
pip install numpy
pip install matplotlib
EOF

echo "anaconda installation complete."
