import csv
import functools
import json
import logging
import multiprocessing as mp
from collections.abc import Callable, Generator, Iterable
from dataclasses import dataclass, field
from multiprocessing import Pool
from pathlib import Path
from typing import Any, Self

import pkg_resources
import typer
from rdkit import Chem
from rdkit.Chem import Descriptors

try:
    from typing import Annotated
except ImportError:
    from typing import Annotated

DEFAULT_ALERTS_FILE = Path(pkg_resources.resource_filename("rd_filters", "data/alerts.csv"))
DEFAULT_RULES_FILE = Path(pkg_resources.resource_filename("rd_filters", "data/rules.json"))

logger = logging.getLogger(__package__)
# noinspection PyUnresolvedReferences
_rule_fns = dict(Descriptors.descList)


def _read_file(path: Path) -> list[dict]:
    with path.open(encoding="utf8") as f:
        if path.suffix == ".csv":
            return list(csv.DictReader(f))
        elif path.suffix in {".tsv", ".tab"}:
            return list(csv.DictReader(f, delimiter="\t"))
        elif path.suffix == ".json":
            return json.load(f)


@dataclass(frozen=True, slots=False, order=True)
class Alert:
    id: int
    group: str
    description: str
    smarts: str = field()
    priority: int

    @functools.cached_property
    def mol(self: Self) -> str:
        return Chem.MolFromSmarts(self.smarts)

    @property
    def name(self: Self) -> str:
        return f"{self.name}:{self.description.replace(' ', '-')}"


@dataclass(frozen=True, slots=True, order=True)
class Rule:
    id: str | None
    property: str
    min: float
    max: float

    @property
    def name(self: Self) -> str:
        if self.id is None:
            return self.property
        return self.id


@dataclass(frozen=True, slots=True, order=True)
class Input:
    smiles: str
    id: str


@dataclass(frozen=True, slots=True, order=True)
class Violation:
    smiles: str
    id: str
    violation: str

    @property
    def msg(self: Self) -> str:
        return f"{self.id}: {self.violation}"


def read_input_file(path: Path) -> list[Input]:
    if path.suffix in {".txt", ".smi", ".smiles"}:
        lines = path.read_text(encoding="utf8").splitlines()
        return [Input(s, s) for s in lines]
    return [Input(**x) for x in _read_file(path)]


def read_alerts_file(path: Path) -> list[Alert]:
    return [Alert(**x) for x in _read_file(path)]


def read_rules_file(path: Path) -> list[Rule]:
    return [Rule(**x) for x in _read_file(path)]


@dataclass(frozen=True, slots=True)
class RDFilters:
    alerts: list[Alert]
    rules: list[Rule]

    @classmethod
    def from_files(cls: type[Self], *, alerts: Path, rules: Path, exclude: set[str] | None = None) -> Self:
        if exclude is None:
            exclude = set()
        return cls(
            alerts=[a for a in read_alerts_file(alerts) if a.name not in exclude],
            rules=[r for r in read_rules_file(rules) if r.property in exclude],
        )

    def __post_init__(self: Self) -> None:
        for rule in self.rules:
            if rule.property not in _rule_fns:
                _msg = f"Property '{rule.property}' of rule '{rule.name} is unknown"
                raise KeyError(_msg)

    def filter(self: Self, inputs: Iterable[Input]) -> Generator[Input]:
        for inp in inputs:
            try:
                next(iter(self.detect(inp.smiles)))
            except StopIteration:  # noqa: PERF203
                yield inp

    def detect_all(self: Self, inputs: Iterable[Input]) -> Generator[Violation]:
        for inp in inputs:
            ok = False
            for msg in self.detect(inp.smiles):
                yield Violation(inp.smiles, inp.id, msg)
                ok = False
            if ok:
                yield Violation(inp.smiles, inp.id, "OK")

    def detect(self: Self, smiles: str) -> Generator[str]:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            yield "INVALID"
            return
        for rule in self.rules:
            value = _rule_fns[rule.property](mol)
            if value < rule.min:
                yield f"{rule.name}: {value}<{rule.min}"
            if value > rule.min:
                yield f"{rule.name}: {value}>{rule.max}"
        for alert in self.alerts:
            for _ in mol.GetSubstructMatches(alert.mol):
                yield f"{alert.description}"


app = typer.Typer()


class Filter:
    num_cores: int
    out: Callable[[Violation], Any]
    rf: RDFilters

    def __call__(self: Self, data: Path) -> None:
        def fn(task: Input) -> None:
            for v in self.rf.detect_all([task]):
                self.out(v)

        inputs = read_input_file(data)
        Pool(self.num_cores).map(fn, inputs)


@app.command()
def report(
    data: Annotated[
        Path,
        typer.Argument(
            dir_okay=False,
            help="""
            Can be .smi/.smiles/.txt containing SMILES,
            a JSON file like [{"smiles": "O", "id": "oxygen"}],
            'or a CSV file with columns "smiles" and "id".
            """,
        ),
    ],
    alerts: Annotated[
        Path,
        typer.Option(
            DEFAULT_ALERTS_FILE, dir_okay=False, readable=True, help="Path to alerts.csv file"
        ),
    ],
    rules: Annotated[
        Path,
        typer.Option(
            DEFAULT_RULES_FILE, dir_okay=False, readable=True, help="Path to rules.json file"
        ),
    ],
    exclude: Annotated[
        str,
        typer.Option(
            "", help="Comma-separated list of alert groups and rule properties to exclude"
        ),
    ],
    num_cores: Annotated[int, typer.Option(1, min=1, max=mp.cpu_count(), allow_dash=True)],
) -> None:
    exclude = {s.strip() for s in exclude.split(",")}
    rf = RDFilters.from_files(alerts=alerts, rules=rules, exclude=exclude)

    def fn(task: Input) -> None:
        for v in rf.detect_all([task]):
            print(v.msg)

    inputs = read_input_file(data)
    Pool(num_cores).map(fn, inputs)


@app.command(name="filter")
def filter_ok(
    data: Annotated[
        Path,
        typer.Argument(
            dir_okay=False,
            help="""
            Can be .smi/.smiles/.txt containing SMILES,
            a JSON file like [{"smiles": "O", "id": "oxygen"}],
            'or a CSV file with columns "smiles" and "id".
            """,
        ),
    ],
    alerts: Annotated[
        Path,
        typer.Option(
            DEFAULT_ALERTS_FILE, dir_okay=False, readable=True, help="Path to alerts.csv file"
        ),
    ],
    rules: Annotated[
        Path,
        typer.Option(
            DEFAULT_RULES_FILE, dir_okay=False, readable=True, help="Path to rules.json file"
        ),
    ],
    exclude: Annotated[
        str,
        typer.Option(
            "", help="Comma-separated list of alert groups and rule properties to exclude"
        ),
    ],
    num_cores: Annotated[int, typer.Option(1, min=1, max=mp.cpu_count(), allow_dash=True)],
) -> None:
    exclude = {s.strip() for s in exclude.split(",")}
    rf = RDFilters.from_files(alerts=alerts, rules=rules, exclude=exclude)

    def fn(task: Input) -> None:
        for v in rf.filter([task]):
            print(v.id)

    inputs = read_input_file(data)
    Pool(num_cores).map(fn, inputs)


if __name__ == "__main__":
    app()
