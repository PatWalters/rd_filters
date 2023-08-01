import csv
import json
import logging
import multiprocessing as mp
import sys
from abc import abstractmethod
from collections.abc import Callable, Generator, Iterable
from dataclasses import asdict, dataclass, field
from functools import cached_property, partialmethod
from multiprocessing import Pool
from pathlib import Path
from typing import Any, Generic, Literal, Self, TypeVar

import pkg_resources
import typer
from rdkit import Chem
from rdkit.Chem import Descriptors

try:
    from typing import Annotated, Unpack
except ImportError:
    from typing import Annotated, Unpack


_logger = logging.getLogger(__package__)
# noinspection PyUnresolvedReferences
_rule_fns = dict(Descriptors.descList)
_empty_set = frozenset()
_T_co = TypeVar("_T_co", covariant=True)


def _read_file(path: Path) -> list[dict]:
    with path.open(encoding="utf8") as f:
        if path.suffix == ".csv":
            return list(csv.DictReader(f))
        elif path.suffix in {".tsv", ".tab"}:
            return list(csv.DictReader(f, delimiter="\t"))
        elif path.suffix == ".json":
            return json.load(f)


@dataclass(frozen=True, slots=False, order=True, kw_only=True)
class Alert:
    name: str
    group: str
    description: str
    smarts: str = field()
    priority: int

    @cached_property
    def mol(self: Self) -> str:
        return Chem.MolFromSmarts(self.smarts)


@dataclass(frozen=True, slots=True, order=True, kw_only=True)
class Rule:
    name: str
    property: str
    min: float
    max: float


@dataclass(frozen=True, slots=False, order=True, kw_only=True)
class Input:
    smiles: str
    id: str

    @cached_property
    def mol(self: Self) -> Chem:
        return Chem.MolFromSmiles(self.smiles)

    def __str__(self: Self) -> str:
        return str(asdict(self))


@dataclass(frozen=True, slots=True, order=True, kw_only=True)
class Violation(Input):

    def __str__(self: Self) -> str:
        return str(asdict(self))

    @classmethod
    def invalid(cls: type[Self], inp: Input) -> Self:
        return cls(id=inp.id, smiles=inp.smiles, type="invalid", name="INVALID SMILES", data={})

    @classmethod
    def alert(cls: type[Self], inp: Input, alert: Alert, **kwargs: Unpack[dict]) -> Self:
        data = asdict(alert) | kwargs
        return cls(id=inp.id, smiles=inp.smiles, type="alert", name=alert.name, data=data)

    @classmethod
    def rule(cls: type[Self], inp: Input, rule: Rule, **kwargs: Unpack[dict]) -> Self:
        data = asdict(rule) | kwargs
        return cls(id=inp.id, smiles=inp.smiles, type="rule", name=rule.name, data=data)

    id: str
    smiles: str
    type: Literal["rule"] | Literal["alert"] | Literal["invalid"]
    name: str
    data: dict[str, Any]


@dataclass(frozen=True, slots=True)
class RDFilters:
    alerts: list[Alert]
    rules: list[Rule]

    def __post_init__(self: Self) -> None:
        for rule in self.rules:
            if rule.property not in _rule_fns:
                _msg = f"Property '{rule.property}' of rule '{rule.name} is unknown"
                raise KeyError(_msg)

    def filter(self: Self, inputs: Iterable[Input]) -> Generator[Input]:
        for inp in inputs:
            try:
                next(iter(self.detect(inp)))
            except StopIteration:  # noqa: PERF203
                yield inp

    def detect_all(self: Self, inputs: Iterable[Input]) -> Generator[Violation]:
        for inp in inputs:
            yield from self.detect(inp)

    def detect(self: Self, inp: Input) -> Generator[Violation]:
        if inp.mol is None:
            yield Violation.invalid(inp)
            return
        for rule in self.rules:
            value = _rule_fns[rule.property](inp.mol)
            if value < rule.min:
                yield Violation.rule(inp, rule, value=value, msg=f"{value} < {rule.min}")
            if value > rule.min:
                yield Violation.rule(inp, rule, value=value, msg=f"{value} > {rule.max}")
        for alert in self.alerts:
            for _ in inp.mol.GetSubstructMatches(alert.mol):
                yield Violation.alert(inp, alert)


class Io:
    DEFAULT_ALERTS_FILE = Path(pkg_resources.resource_filename("rd_filters", "data/alerts.csv"))
    DEFAULT_RULES_FILE = Path(pkg_resources.resource_filename("rd_filters", "data/rules.json"))

    @staticmethod
    def read_input_file(path: Path) -> list[Input]:
        if path.suffix in {".txt", ".smi", ".smiles"}:
            lines = path.read_text(encoding="utf8").splitlines()
            return [Input(s, s) for s in lines]
        return [Input(**x) for x in _read_file(path)]


    @staticmethod
    def read_alerts_file(path: Path) -> list[Alert]:
        return [
            Alert(**({"name": x["description"]} | x))
            for x in _read_file(path)
        ]


    @staticmethod
    def read_rules_file(path: Path) -> list[Rule]:
        return [
            Rule(**({"name": x["property"]} | x))
            for x in _read_file(path)
        ]


@dataclass(frozen=True, kw_only=True)
class Processor(Generic[_T_co]):

    @abstractmethod
    def task(self: Self, rf: RDFilters, smiles: Input) -> _T_co:
        raise NotImplementedError()

    alerts: Path = Io.DEFAULT_ALERTS_FILE
    rules: Path = Io.DEFAULT_RULES_FILE
    exclude: set[str] | frozenset[str] = _empty_set
    num_cores: int = 1

    @cached_property
    def rf(self: Self) -> RDFilters:
        return RDFilters(
            alerts=[a for a in Io.read_alerts_file(self.alerts) if a.name not in self.exclude],
            rules=[r for r in Io.read_rules_file(self.rules) if r.property in self.exclude],
        )

    def __call__(self: Self, data: Path) -> None:
        inputs = Io.read_input_file(data)
        Pool(self.num_cores).map(partialmethod(self.task, rf=self.rf), inputs)


def _violation_msg_writer(v: Violation) -> None:
    if msg := v.data.get("msg"):
        sys.stdout.write(f"{v.id} : {v.name} : {msg}")
    else:
        sys.stdout.write(f"{v.id} : {v.name}")

def _violation_json_writer(v: Violation) -> None:
    sys.stdout.write(str(v))

def _input_id_writer(inp: Input) -> None:
    sys.stdout.write(inp.id)

def _input_json_writer(inp: Input) -> None:
    sys.stdout.write(str(inp))


@dataclass(frozen=True, kw_only=True)
class ViolationPrinter(Processor):
    writer: Callable[[Violation], Any] = _violation_msg_writer

    def task(self: Self, rf: RDFilters, input_: Input) -> None:
        for v in rf.detect_all([input_]):
            self.writer(v)


@dataclass(frozen=True, kw_only=True)
class FilterPrinter(Processor):
    writer: Callable[[Input], Any] = _input_id_writer

    def task(self: Self, rf: RDFilters, input_: Input) -> None:
        for v in rf.filter([input_]):
            self.writer(v)


class _Args:
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
    ] = Path("input.smiles")
    alerts: Annotated[
        Path,
        typer.Option(
            Io.DEFAULT_ALERTS_FILE, dir_okay=False, readable=True, help="Path to alerts.csv file"
        ),
    ] = Io.DEFAULT_ALERTS_FILE
    rules: Annotated[
        Path,
        typer.Option(
            Io.DEFAULT_RULES_FILE, dir_okay=False, readable=True, help="Path to rules.json file"
        ),
    ] = Io.DEFAULT_RULES_FILE
    exclude: Annotated[
        str,
        typer.Option(
            "", help="Comma-separated list of alert groups and rule properties to exclude"
        ),
    ] = ""
    violation_as_json: Annotated[
        bool,
        typer.Option(
            False, "--json", help="Write full JSON data per violation instead of message"
        )
    ] = False
    input_as_json: Annotated[
        bool,
        typer.Option(
            False, "--json", help="Write JSON string that includes SMILES, rather than just the ID"
        )
    ] = False
    num_cores: Annotated[int, typer.Option(1, min=1, max=mp.cpu_count(), allow_dash=True)] = 1


app = typer.Typer()

@app.command(
    name="report",
    short_help="",
    help="",
)
def report(
    data = _Args.data,  # noqa: ANN001
    alerts = _Args.alerts,  # noqa: ANN001
    rules = _Args.rules,  # noqa: ANN001
    exclude = _Args.exclude,  # noqa: ANN001
    as_json = _Args.violation_as_json,  # noqa: ANN001
    num_cores = _Args.num_cores,  # noqa: ANN001
) -> None:
    processor = ViolationPrinter(
        alerts=alerts,
        rules=rules,
        exclude={s.strip() for s in exclude.split(",")},
        num_cores=num_cores,
        writer=_violation_json_writer if as_json else _violation_msg_writer,
    )
    processor(data)

@app.command(
    name="filter",
    short_help="",
    help="",
)
def filter_ok(
    data = _Args.data,  # noqa: ANN001
    alerts = _Args.alerts,  # noqa: ANN001
    rules = _Args.rules,  # noqa: ANN001
    exclude = _Args.exclude,  # noqa: ANN001
    as_json = _Args.input_as_json,  # noqa: ANN001
    num_cores = _Args.num_cores,  # noqa: ANN001
) -> None:
    processor = FilterPrinter(
        alerts=alerts,
        rules=rules,
        exclude={s.strip() for s in exclude.split(",")},
        num_cores=num_cores,
        writer=_input_json_writer if as_json else _input_id_writer
    )
    processor(data)


if __name__ == "__main__":
    app()
