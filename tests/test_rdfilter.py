import pytest

from rd_filters import rd_filters

# these just stop on the first filter
_lint_data = {
    "C1N=C1": {"aziridine-like N in 3-membered ring", "aziridine_diazirine"},
    "NN=N": {
        "Azo",
        "Oxygen-nitrogen single bond",
        "acyclic N-,=N and not N bound to carbonyl or sulfone",
        "azo_A(324)",
        "azo_amino",
        "diazo group",
        "hydrazine",
    },
    "CN=C": {"imine"},  # WAS "OK"
    "CN=CC": {r"Filter39_imine", "imine", "acyclic imines"},  # WAS "OK"
    "CN=NC": {r"Filter5_azo"},
    "CN=N": {r"Filter5_azo"},
    "N=N": {r"Filter5_azo"},
    "{Fe}C": {r"Filter9_metal"},
    "{FeH}": {r"Filter9_metal"},
    "{Fe}": {"OK"},
    "CN=O": {r"Filter12_nitroso"},
    "N=O": {"Filter12_nitroso"},
    "S=P": {"Filter13_PS_double_bond"},
    "S=PC": {"Filter13_PS_double_bond"},
    "CC=OSC": {"Filter29_thioester"},
    "CC=SSC": {"Filter29_thioester"},
    "CC=NSC": {"Filter29_thioester"},
    "CC=OS": {"Filter29_thioester"},
    "CC=SS": {"Filter29_thioester"},
    "CC=NS": {"Filter29_thioester"},
    "CC=OCCBr": {"Filter26_alkyl_halide"},
    "CC=OCCI": {"Filter26_alkyl_halide"},
    "C=OCCI": {"Filter26_alkyl_halide"},
    "NC=OCCI": {"Filter26_alkyl_halide"},
    "CC=NOS=ON": {"Filter18_oxime_ester"},
    "NC=NOP=ON": {"Filter89_hydroxylamine"},
    "C=NOP=ON": {"Filter18_oxime_ester"},
    "SO": {"Filter31_so_bond"},
    "SOC": {"Filter31_so_bond"},
    "c1csocc1": {"OK"},
    "OO": {"Filter32_oo_bond"},
    "COO": {"Filter32_oo_bond"},
    "c1coocc1": {"OK"},
    "CC=OC=C": {"Filter44_michael_acceptor2"},
    "C=OC=C": {"Filter38_aldehyde"},
    "OC=OC=C": {"Filter44_michael_acceptor2"},
    "C1=OC=CC=OC=C1": {"Filter53_para_quinones"},
    "CIC": {"Filter49_halogen"},
    "CICC": {"Filter49_halogen"},
    "NC=NO": {"Filter89_hydroxylamine"},
    "C=ONOS": {"Filter31_so_bond"},
}


class Test:
    def test_detect(self):
        alerts = rd_filters.read_alerts_file(rd_filters.DEFAULT_ALERTS_FILE)
        rf = rd_filters.RDFilters(alerts, [])
        for smiles, expected in _lint_data.items():
            actual = set(rf.detect(smiles))
            assert actual == expected, f"Invalid result for '{smiles}'"
            # assert {r.msg for r in actual} == [f"{y}: {x}" for x, y in _inpharmatica_data]


if __name__ == "__main__":
    pytest.main()
