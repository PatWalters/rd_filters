import csv
import functools
import json
import logging
import multiprocessing as mp
from collections.abc import Generator, Iterable
from dataclasses import dataclass, field
from multiprocessing import Pool
from pathlib import Path
from typing import Optional

import pkg_resources
import typer
from rdkit import Chem
from rdkit.Chem.Descriptors import TPSA, MolLogP, MolWt, NumHAcceptors, NumHDonors
from rdkit.Chem.rdMolDescriptors import CalcNumRotatableBonds

logger = logging.getLogger(__package__)
_alerts_file = pkg_resources.resource_filename("rd_filters", "data/alerts.csv")
_rules_file = pkg_resources.resource_filename("rd_filters", "data/rules.json")
_rule_fns = {
    "MW": MolWt,
    "LogP": MolLogP,
    "HBD": NumHDonors,
    "HBA": NumHAcceptors,
    "TPSA": TPSA,
    "Rot": CalcNumRotatableBonds,
}


def _read_file(path: Path) -> list[dict]:
    with path.open(encoding="utf8") as f:
        if path.suffix == ".csv":
            return list(csv.DictReader(f))
        elif path.suffix in {".tsv", ".tab"}:
            return list(csv.DictReader(f, delimiter="\t"))
        elif path.suffix == ".json":
            return json.load(f)


@dataclass(frozen=True, repr=True, order=True)
class Alert:
    id: int
    group: str
    description: str
    smarts: str = field()
    priority: int

    @functools.cached_property
    def chem(self) -> str:
        return Chem.MolFromSmarts(self.smarts)

    @property
    def name(self) -> str:
        return f"{self.name}:{self.description.replace(' ', '-')}"


@dataclass(frozen=True, repr=True, order=True)
class Rule:
    id: str | None
    property: int
    min: float
    max: float

    @property
    def name(self):
        if self.id is None:
            return self.property
        else:
            return self.id


@dataclass(frozen=True, repr=True, order=True)
class Input:
    smiles: str
    id: str


@dataclass(frozen=True, repr=True, order=True)
class Violation:
    smiles: str
    id: str
    violation: str

    @property
    def msg(self) -> str:
        return f"{self.id}: {self.violation}"


def _read_input_file(path: Path) -> list[Input]:
    if path.suffix in {".txt", ".smi", ".smiles"}:
        lines = path.read_text(encoding="utf8").splitlines()
        return [Input(s, s) for s in lines]
    return [Input(**x) for x in _read_file(path)]


def _read_alerts_file(path: Path) -> list[Alert]:
    return [Alert(**x) for x in _read_file(path)]


def _read_rules_file(path: Path) -> list[Rule]:
    return [Rule(**x) for x in _read_file(path)]


@dataclass(frozen=True)
class RDFilters:
    alerts: list[Alert]
    rules: list[Rule]

    def __post_init__(self):
        for rule in self.rules:
            if rule.property not in _rule_fns:
                raise KeyError(f"Property '{rule.property}' of rule '{rule.name} is unknown")

    def detect_all(self, inputs: Iterable[Input]) -> Generator[Violation, None, None]:
        for inp in inputs:
            for msg in self.detect(inp.smiles):
                yield Violation(inp.smiles, inp.id, msg)
            else:
                yield Violation(inp.smiles, inp.id, "OK")

    def detect(self, smiles: str) -> Generator[str, None, None]:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            yield "invalid"
            return
        for rule in self.rules:
            value = _rule_fns[rule.property](mol)
            if value < rule.min:
                yield f"{rule.name}: {value}<{rule.min}"
            if value > rule.min:
                yield f"{rule.name}: {value}>{rule.max}"
        for alert in self.alerts:
            for match in mol.GetSubstructMatches(alert.smarts):
                yield f"{alert.description}: '{match}'"


app = typer.Typer()


@app.command()
def report(
    data: Path = typer.Argument(
        dir_okay=False,
        help="""
        Can be .smi/.smiles/.txt containing SMILES,
        a JSON file like [{"smiles": "O", "id": "oxygen"}],
        'or a CSV file with columns "smiles" and "id".
        """,
    ),
    alerts: Path = typer.Option(_alerts_file, dir_okay=False, help="Path to alerts.csv file"),
    rules: Path = typer.Option(_rules_file, dir_okay=False, help="Path to rules.json file"),
    exclude: str = typer.Option(
        "", help="Comma-separated list of alert groups and rule properties to exclude"
    ),
    num_cores: int = typer.Option(1, min=1, max=mp.cpu_count()),
):
    exclude = {s.strip() for s in exclude.split(",")}
    alerts = [a for a in _read_alerts_file(alerts) if a.name not in exclude]
    rules = [r for r in _read_rules_file(rules) if r.property in exclude]
    rf = RDFilters(alerts, rules)

    def fn(task):
        for v in rf.detect_all([task]):
            print(v.msg)

    inputs = _read_input_file(data)
    Pool(num_cores).map(fn, inputs)


@app.command(name="filter")
def filter_ok(
    data: Path = typer.Argument(
        dir_okay=False,
        help="""
        Can be .smi/.smiles/.txt containing SMILES,
        a JSON file like [{"smiles": "O", "id": "oxygen"}],
        'or a CSV file with columns "smiles" and "id".
        """,
    ),
    alerts: Path = typer.Option(_alerts_file, dir_okay=False, help="Path to alerts.csv file"),
    rules: Path = typer.Option(_rules_file, dir_okay=False, help="Path to rules.json file"),
    exclude: str = typer.Option(
        "", help="Comma-separated list of alert groups and rule properties to exclude"
    ),
    num_cores: int = typer.Option(1, min=1, max=mp.cpu_count()),
):
    exclude = {s.strip() for s in exclude.split(",")}
    alerts = [a for a in _read_alerts_file(alerts) if a.name not in exclude]
    rules = [r for r in _read_rules_file(rules) if r.property in exclude]
    rf = RDFilters(alerts, rules)

    def fn(task):
        for v in rf.detect_all([task]):
            if v.violation == "OK":
                print(v.id)
                break

    inputs = _read_input_file(data)
    Pool(num_cores).map(fn, inputs)


if __name__ == "__main__":
    app()
