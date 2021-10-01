#!/usr/bin/env python3

import sys
from rdkit import Chem
from rdkit.Chem.Descriptors import MolWt, MolLogP, NumHDonors, NumHAcceptors, TPSA
from rdkit.Chem.rdMolDescriptors import CalcNumRotatableBonds
import multiprocessing as mp
from multiprocessing import Pool
import time
import pandas as pd
import os
import json
from docopt import docopt
import pkg_resources

cmd_str = """Usage:
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


def read_rules(rules_file_name):
    """
    Read rules from a JSON file
    :param rules_file_name: JSON file name
    :return: dictionary corresponding to the contents of the JSON file
    """
    with open(rules_file_name) as json_file:
        try:
            rules_dict = json.load(json_file)
            return rules_dict
        except json.JSONDecodeError:
            print(f"Error parsing JSON file {rules_file_name}")
            sys.exit(1)


def write_rules(rule_dict, file_name):
    """
    Write configuration to a JSON file
    :param rule_dict: dictionary with rules
    :param file_name: JSON file name
    :return: None
    """
    ofs = open(file_name, "w")
    ofs.write(json.dumps(rule_dict, indent=4, sort_keys=True))
    print(f"Wrote rules to {file_name}")
    ofs.close()


def default_rule_template(alert_list, file_name):
    """
    Build a default rules template
    :param alert_list: list of alert set names
    :param file_name: output file name
    :return: None
    """
    default_rule_dict = {
        "MW": [0, 500],
        "LogP": [-5, 5],
        "HBD": [0, 5],
        "HBA": [0, 10],
        "TPSA": [0, 200],
        "Rot": [0, 10]
    }
    for rule_name in alert_list:
        if rule_name == "Inpharmatica":
            default_rule_dict["Rule_" + rule_name] = True
        else:
            default_rule_dict["Rule_" + rule_name] = False
    write_rules(default_rule_dict, file_name)


def get_config_file(file_name, environment_variable):
    """
    Read a configuration file, first look for the file, if you can't find
    it there, look in the directory pointed to by environment_variable
    :param file_name: the configuration file
    :param environment_variable: the environment variable
    :return: the file name or file_path if it exists otherwise exit
    """
    if os.path.exists(file_name):
        return file_name
    else:
        config_dir = os.environ.get(environment_variable)
        if config_dir:
            config_file_path = os.path.join(os.path.sep, config_dir, file_name)
            if os.path.exists(config_file_path):
                return config_file_path

    error_list = [f"Could not file {file_name}"]
    if config_dir:
        err_str = f"Could not find {config_file_path} based on the {environment_variable}" + \
                  "environment variable"
        error_list.append(err_str)
    error_list.append(f"Please check {file_name} exists")
    error_list.append(f"Or in the directory pointed to by the {environment_variable} environment variable")
    print("\n".join(error_list))
    sys.exit(1)


class RDFilters:
    def __init__(self, rules_file_name):
        good_name = get_config_file(rules_file_name, "FILTER_RULES_DIR")
        self.rule_df = pd.read_csv(good_name)
        # make sure there wasn't a blank line introduced
        self.rule_df = self.rule_df.dropna()
        self.rule_list = []

    def build_rule_list(self, alert_name_list):
        """
        Read the alerts csv file and select the rule sets defined in alert_name_list
        :param alert_name_list: list of alert sets to use
        :return:
        """
        self.rule_df = self.rule_df[self.rule_df.rule_set_name.isin(alert_name_list)]
        tmp_rule_list = self.rule_df[["rule_id", "smarts", "max", "description"]].values.tolist()
        for rule_id, smarts, max_val, desc in tmp_rule_list:
            smarts_mol = Chem.MolFromSmarts(smarts)
            if smarts_mol:
                self.rule_list.append([smarts_mol, max_val, desc])
            else:
                print(f"Error parsing SMARTS for rule {rule_id}", file=sys.stderr)

    def get_alert_sets(self):
        """
        :return: a list of unique rule set names
        """
        return self.rule_df.rule_set_name.unique()

    def evaluate(self, lst_in):
        """
        Evaluate structure alerts on a list of SMILES
        :param lst_in: input list of [SMILES, Name]
        :return: list of alerts matched or "OK"
        """
        smiles, name = lst_in
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return [smiles, name, 'INVALID', -999, -999, -999, -999, -999]
        desc_list = [MolWt(mol), MolLogP(mol), NumHDonors(mol), NumHAcceptors(mol), TPSA(mol),
                     CalcNumRotatableBonds(mol)]
        for row in self.rule_list:
            patt, max_val, desc = row
            if len(mol.GetSubstructMatches(patt)) > max_val:
                return [smiles, name] + [desc + " > %d" % (max_val)] + desc_list
        return [smiles, name] + ["OK"] + desc_list


def main():
    cmd_input = docopt(cmd_str)
    alert_file_name = cmd_input.get("--alerts") or pkg_resources.resource_filename('rd_filters',
                                                                                   "data/alert_collection.csv")
    rf = RDFilters(alert_file_name)

    if cmd_input.get("template"):
        template_output_file = cmd_input.get("--out")
        default_rule_template(rf.get_alert_sets(), template_output_file)

    elif cmd_input.get("filter"):
        input_file_name = cmd_input.get("--in")
        rules_file_name = cmd_input.get("--rules") or pkg_resources.resource_filename('rd_filters', "data/rules.json")
        rules_file_path = get_config_file(rules_file_name, "FILTER_RULES_DATA")
        prefix_name = cmd_input.get("--prefix")
        num_cores = cmd_input.get("--np") or mp.cpu_count()
        num_cores = int(num_cores)

        print("using %d cores" % num_cores, file=sys.stderr)
        start_time = time.time()
        p = Pool(num_cores)
        input_data = [x.split() for x in open(input_file_name)]
        input_data = [x for x in input_data if len(x) == 2]
        rule_dict = read_rules(rules_file_path)

        rule_list = [x.replace("Rule_", "") for x in rule_dict.keys() if x.startswith("Rule") and rule_dict[x]]
        rule_str = " and ".join(rule_list)
        print(f"Using alerts from {rule_str}", file=sys.stderr)
        rf.build_rule_list(rule_list)
        res = list(p.map(rf.evaluate, input_data))
        df = pd.DataFrame(res, columns=["SMILES", "NAME", "FILTER", "MW", "LogP", "HBD", "HBA", "TPSA", "Rot"])
        df_ok = df[
            (df.FILTER == "OK") &
            df.MW.between(*rule_dict["MW"]) &
            df.LogP.between(*rule_dict["LogP"]) &
            df.HBD.between(*rule_dict["HBD"]) &
            df.HBA.between(*rule_dict["HBA"]) &
            df.TPSA.between(*rule_dict["TPSA"]) &
            df.Rot.between(*rule_dict["Rot"])
            ]
        output_smiles_file = prefix_name + ".smi"
        output_csv_file = prefix_name + ".csv"
        df_ok[["SMILES", "NAME"]].to_csv(f"{output_smiles_file}", sep=" ", index=False, header=False)
        print(f"Wrote SMILES for molecules passing filters to {output_smiles_file}", file=sys.stderr)
        df.to_csv(f"{prefix_name}.csv", index=False)
        print(f"Wrote detailed data to {output_csv_file}", file=sys.stderr)

        num_input_rows = df.shape[0]
        num_output_rows = df_ok.shape[0]
        fraction_passed = "%.1f" % (num_output_rows / num_input_rows * 100.0)
        print(f"{num_output_rows} of {num_input_rows} passed filters {fraction_passed}%", file=sys.stderr)
        elapsed_time = "%.2f" % (time.time() - start_time)
        print(f"Elapsed time {elapsed_time} seconds", file=sys.stderr)


if __name__ == "__main__":
    main()
