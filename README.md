# rd_filters

This script provides a simple means of applying the
functional group alerts from the ChEMBL database and property alerts from the RDKit
to a set of compounds. As of ChEMBL 23, the database table **structural_alerts**
contains 8 sets of alerts. ChEMBL doesn't have much documentation on the different alert sets.

| Rule Set                                                | Number of Alerts |
| ------------------------------------------------------- | ---------------: |
| BMS                                                     |              180 |
| Dundee                                                  |              105 |
| Glaxo                                                   |               55 |
| Inpharmatica                                            |               91 |
| LINT                                                    |               57 |
| MLSMR                                                   |              116 |
| [PAINS](https://pubs.acs.org/doi/abs/10.1021/jm901137j) |              479 |
| SureChEMBL                                              |              166 |

The SMARTS for some alerts are incompatible with RDKit,
so I edited them. A complete list of the changes I made is in the file **Notes.txt**.

## Prerequisite

- Python 3.7+
- [rdkit](https://www.rdkit.org/docs/Install.html)

## Installation

Run `pip install rd-filters`

Directly from GitHub:
`pip install git+https://github.com/PatWalters/rd_filters.git`

Install locally:

```bash
git clone https://github.com/PatWalters/rd_filters
cd rd_filters
pip install .
```

## Usage

### Examples

For command-line help, run:

```shell
`rd-filter --help`
```

Usage is simple. To filter a file `test.smi` containing line-by-line SMILES:

```bash
rd-filters filter test.smi > filtered.smi
```

If I want to understand the violations that occurred:

```bash
rd-filters report test.smi > report.txt
```

### Config files

2 config files are used. Both can be passed at the command-line.

- `alerts.csv` - the set of structural alerts
  You shouldn't have to mess with this unless you want to add your own structural alerts.
- `rules.json` - the set of other violations.

## Contributing

As always please let me know if you have questions, comments, etc.

Copyright Pat Walters, August 2018
