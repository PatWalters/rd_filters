import pkg_resources
import pytest

from rd_filters import rd_filters

# these just stop on the first filter
_lint_data = [
    ("C1N=C1", "aziridine-like N in 3-membered ring > 0"),
    ("NN=N", "acyclic N-,=N and not N bound to carbonyl or sulfone > 0"),
]

_inpharmatica_data = [
    ("CN=C", "OK"),
    ("CN=CC", "Filter39_imine > 0"),
    ("CN=C(C)", "Filter39_imine > 0"),
    ("CN=NC", "Filter5_azo > 0"),
    ("CN=N", "Filter5_azo > 0"),
    ("N=N", "Filter5_azo > 0"),
    ("[Fe]C", "Filter9_metal > 0"),
    ("[FeH]", "Filter9_metal > 0"),
    ("[Fe]", "OK"),
    ("CN=O", "Filter12_nitroso > 0"),
    ("N=O", "Filter12_nitroso > 0"),
    ("S=P", "Filter13_PS_double_bond > 0"),
    ("S=PC", "Filter13_PS_double_bond > 0"),
    ("CC(=O)SC", "Filter29_thioester > 0"),
    ("CC(=S)SC", "Filter29_thioester > 0"),
    ("CC(=N)SC", "Filter29_thioester > 0"),
    ("CC(=O)S", "Filter29_thioester > 0"),
    ("CC(=S)S", "Filter29_thioester > 0"),
    ("CC(=N)S", "Filter29_thioester > 0"),
    ("CC(=O)CCBr", "Filter26_alkyl_halide > 0"),
    ("CC(=O)CCI", "Filter26_alkyl_halide > 0"),
    ("C(=O)CCI", "Filter26_alkyl_halide > 0"),
    ("NC(=O)CCI", "Filter26_alkyl_halide > 0"),
    ("CC=NOS(=O)N", "Filter18_oxime_ester > 0"),
    ("NC=NOP(=O)N", "Filter89_hydroxylamine > 0"),
    ("C=NOP(=O)N", "Filter18_oxime_ester > 0"),
    ("SO", "Filter31_so_bond > 0"),
    ("SOC", "Filter31_so_bond > 0"),
    ("c1csocc1", "OK"),
    ("OO", "Filter32_oo_bond > 0"),
    ("COO", "Filter32_oo_bond > 0"),
    ("c1coocc1", "OK"),
    ("CC(=O)C=C", "Filter44_michael_acceptor2 > 0"),
    ("C(=O)C=C", "Filter38_aldehyde > 0"),
    ("OC(=O)C=C", "Filter44_michael_acceptor2 > 0"),
    ("C1(=O)C=CC(=O)C=C1", "Filter53_para_quinones > 0"),
    ("CIC", "Filter49_halogen > 0"),
    ("CI(C)C", "Filter49_halogen > 0"),
    ("CC=NOS(=O)N", "Filter18_oxime_ester > 0"),
    ("NC=NO", "Filter89_hydroxylamine > 0"),
    ("C(=O)NOS", "Filter31_so_bond > 0"),
]


def test_lint():
    alerts = rd_filters._read_alerts_file(rd_filters._alerts_file)
    rf = rd_filters.RDFilters(alerts, [])
    for input_, expected in _lint_data:
        actual = list(rf.detect(input_))
        assert actual == [expected]


def test_inpharmatica():
    alerts = rd_filters._read_alerts_file(rd_filters._alerts_file)
    rf = rd_filters.RDFilters(alerts, [])
    input_data = [rd_filters.Input(x, x) for x, y in _inpharmatica_data]
    results = list(rf.detect_all(input_data))
    assert [r.violation for r in results] == [x for x, _ in _inpharmatica_data]


if __name__ == "__main__":
    pytest.main()
