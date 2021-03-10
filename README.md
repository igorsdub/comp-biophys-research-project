Computational Biophysics Research Project
==============================

A descent project structure for doing and sharing science work in Computational Biophysics.

Project Organization
------------

    ├── Makefile           <- Makefile with commands like `make data` or `make train`
    ├── README.md          <- The top-level README for developers using this project.
    ├── config.yaml        <- A configuration file for the scripts.
    |
    ├── data
    │   ├── external       <- Data from third party sources.
    │   ├── interim        <- Intermediate data that has been transformed.
    │   ├── processed      <- The final, canonical data sets for modeling.
    │   └── raw            <- The original, immutable data dump.
    |
    ├── docs               <- A default Sphinx project; see sphinx-doc.org for details
    |
    ├── ccenv.yml          <- A Conda environmet file.
    |
    ├── misc
    |
    ├── notebooks          <- Jupyter notebooks. Naming convention is a number (for ordering),
    │                         the creator's initials, and a short `-` delimited description, e.g.
    │                         `1.0-jqp-initial-data-exploration`.
    |
    ├── output
    |
    ├── pdb
    │   ├── processed
    │   └── raw
    |
    ├── references         <- Data dictionaries, manuals, and all other explanatory materials.
    │
    ├── reports            <- Generated analysis as HTML, PDF, LaTeX, etc.
    │   └── figures        <- Generated graphics and figures to be used in reporting
    |
    ├── src
    │   ├── data           <- Scripts to download or generate data
    │   │   ├── make_dataset.py
    │   │   └── make_pdb.py
    │   └── visualization  <- Scripts to create exploratory and results oriented visualizations
    │       └── visualize.py
    |
    ├── tmp                <- for temporary files
    |
    ├── tox.ini            <- tox file with settings for running tox; see tox.readthedocs.io
    |
    └── .env               <- Stores environment variables for Dotenv module. Do not track with version control.

--------

<p><small>Project based on the <a target="_blank" href="https://drivendata.github.io/cookiecutter-data-science/">cookiecutter data science project template</a>. #cookiecutterdatascience</small></p>
