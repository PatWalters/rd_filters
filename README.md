# rd_filters

This script provides a simple means of applying the 
functional group filters from the ChEMBL database, as well as a number of 
property  filters from the RDKit, to a set 
of compounds.  As of ChEMBL 23,  the database table **structural_alerts**
contains 8 sets of alerts.  ChEMBL doesn't apper to have much in the way of 
documentation on the different alert sets. 

|Rule Set |Number of Alerts|
|---|---:|
|BMS|	180|
|Dundee|	105|
|Glaxo|	55|
|Inpharmatica|	91|
|LINT|	57|
|MLSMR|	116|
|[PAINS](https://pubs.acs.org/doi/abs/10.1021/jm901137j)|	479|
|SureChEMBL|	166|

The SMARTS patterns in a number of these alerts were not compatible with the RDKit
so I edited them.  A complete list of the changes I made is in the file **Notes.txt**. 

## Prerequisite

* At least Python 3.6
* The RDKit, you can find installation instructions [here](https://www.rdkit.org/docs/Install.html).  I'd recommend the conda route.

## Installation

### Directly install from github

`pip install git+https://github.com/PatWalters/rd_filters.git`

### Local install

``` shell
git clone https://github.com/PatWalters/rd_filters
cd rd_filters
pip install .
```

## Usage

The script needs 2 files to operate.

* alert_collection.csv - the set of structural alerts
* rules.json - the configuration file

The script uses the following logic to find alert_collection.csv and rules.json.

1. Use locations specified by the "--alert" (for alerts.csv) and "--rules" (for rules.json) command line arguments.
2. Look in the current directory.
3. Look in the directory pointed to by the FILTER_RULES_DATA environment variable.

I'll provide some examples below to illustrate.  

That's it, at this point you should be good to go. 

### Configuration files

The file **alert_collection.csv** contains alerts.  You shouldn't have to mess with this unless you
want to add your own structural alerts. I think the format is pretty obvious.

The file **rules.json** controls which filters and alerts are used.  You can use the command
below to generate a **rules.json** with the default settings. 

`rd_filters template --out rules.json`

The **rules.json** file looks like this. The values for the properties are the maximum and minimum 
allowed (inclusive).  To set which structural alerts are used, set **true** and **false**. You can 
use multiple alert sets. 
Just edit the file with your [favorite text editor](https://www.gnu.org/software/emacs/).

```javascript
{
    "HBA": [
        0,
        10
    ],
    "HBD": [
        0,
        5
    ],
    "LogP": [
        -5,
        5
    ],
    "MW": [
        0,
        500
    ],
    "Rule_BMS": false,
    "Rule_Dundee": false,
    "Rule_Glaxo": false,
    "Rule_Inpharmatica": true,
    "Rule_LINT": false,
    "Rule_MLSMR": false,
    "Rule_PAINS": false,
    "Rule_SureChEMBL": false,
    "TPSA": [
        0,
        200
    ]
}
```

#### Examples

First off, you're going to want to copy **alert_collection.csv** and 
**rules.json** to a directory and set the FILTER_RULES_DATA environment
variable to point to that directory.  If you are using a bash-ish shell
and the files are in /home/elvis/data that would be:

`export FILTER_RULES_DATA=/home/elvis/data`

If you type

`rd_filters -h`

you'll see this:

``` shell
Usage:
rd_filters filter --in INPUT_FILE --prefix PREFIX [--rules RULES_FILE_NAME] [--alerts ALERT_FILE_NAME][--np NUM_CORES]
rd_filters template --out TEMPLATE_FILE [--rules RULES_FILE_NAME]

Options:
--in INPUT_FILE input file name
--prefix PREFIX prefix for output file names
--rules RULES_FILE_NAME name of the rules JSON file
--alerts ALERTS_FILE_NAME name of the structural alerts file
--np NUM_CORES the number of cpu cores to use (default is all)
--out TEMPLATE_FILE parameter template file name
"""
```

The basic operation is pretty simple. If I want to filter a file called test.smi and 
I want my output files to start with "out", I could do something like this:

`rd_filters filter --in test.smi --prefix out`

This will create 2 files
* **out.smi** - contains the SMILES strings and molecule names for all of the compounds
passing the filters
* **out.csv** - contains calculated property values and a listing alerts triggered
by a molecule

By default, this script runs in parallel and uses all available processors.  To
change this value, use the --np flag. 

`rd_filters filter --in test.smi --prefix out --np 4`

As mentioned above, alternate rules files or alerts files can be specified on
the command line.

``` shell
rd_filters filter --in test.smi --prefix out --rules myrules.json
rd_filters filter --in test.smi --prefix out --alerts myalerts.csv
rd_filters filter --in test.smi --prefix out --rules myrules.json --alerts myalerts.csv
```

A new default rules template file can be generated using the **template** option.

`rdfilters.py template --out myrules.json`

As always please let me know if you have questions, comments, etc. 

Pat Walters, August 2018
