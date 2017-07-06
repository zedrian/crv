CRV
===
Chromatography results verifier
___

## Requirements
* **Python 3** with installed `pandas` package
* **CFM-ID** executables (`cfm-id.exe` and others) that should be placed in `cfm-id` folder

## TODO
* it's necessary to show SMILE for all found fragments
* provide mechanism of chemical elements' selection for mass check
* provide mechanism for database generation (using SMILEs and `cfm-predict`)
* provide full-power work with database: initialization, add / remove elements, data analysis (searching for common and unique fragments)
* implement per-compound analysis (from sample to sample)

## N.B.
* empty lines in `all-compounds.csv` are possible and it's okay
* script should be able to work with both folder and folder-with-folders
* `cfm-annotate` can be used for generating SMILEs for fragments (in addition to `cfm-predict`)
* `cfm-id` should be more effective than `cfm-id-precomputed`, but needs correct training model with `cfm-train`