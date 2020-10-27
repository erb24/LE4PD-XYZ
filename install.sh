#!/bin/bash

BD=$PWD

pip install -e .
cd LE4PD/util/
rm -rfv __pycache__
rm -rfv *.so

cd $BD
