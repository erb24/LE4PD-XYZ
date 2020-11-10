from setuptools import setup

setup(
    name = 'LE4PD',
    version = '0.1',
    author = 'Eric Beyerle',
    description = 'Python API for Protein Dynamics using the Langevin Formalism (LE4PD)',
    url = 'https://github.com/erb24/LE4PD',

    packages = ['LE4PD'],
    install_requires=[
        'numpy',
        'scipy',
        'matplotlib',
        'mdtraj',
	'physt'
    ],
    zip_safe = False
)
