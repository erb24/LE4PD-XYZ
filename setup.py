from setuptools import setup

setup(
    name = 'LE4PD-XYZ',
    version = '0.1',
    author = 'Eric Beyerle',
    description = 'Python API for Protein Dynamics using the Langevin Formalism (LE4PD)',
    url = 'https://github.com/erb24/LE4PD-PCA',

    packages = ['LE4PD-XYZ'],
    install_requires=[
        'numpy',
        'scipy',
        'matplotlib',
	'physt'
    ],
    zip_safe = False
)
