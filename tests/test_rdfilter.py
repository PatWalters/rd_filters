from rd_filters import rd_filters
import pkg_resources

# these just stop on the first filter
test_lint = [
    ('C1N=C1', 'aziridine-like N in 3-membered ring > 0'),
    ('NN=N', 'acyclic N-,=N and not N bound to carbonyl or sulfone > 0'),
    ]

test_inpharmatica = [
    ('CN=C', 'OK'),
    ('CN=CC', 'Filter39_imine > 0'),
    ('CN=C(C)', 'Filter39_imine > 0'),
    ('CN=NC', 'Filter5_azo > 0'),
    ('CN=N', 'Filter5_azo > 0'),
    ('N=N', 'Filter5_azo > 0'),
    ('[Fe]C', 'Filter9_metal > 0'),
    ('[FeH]', 'Filter9_metal > 0'),
    ('[Fe]', 'OK'),
    ('CN=O', 'Filter12_nitroso > 0'),
    ('N=O', 'Filter12_nitroso > 0'),
    ('S=P', 'Filter13_PS_double_bond > 0'),
    ('S=PC', 'Filter13_PS_double_bond > 0'),
    ('CC(=O)SC', 'Filter29_thioester > 0'),
    ('CC(=S)SC', 'Filter29_thioester > 0'),
    ('CC(=N)SC', 'Filter29_thioester > 0'),
    ('CC(=O)S', 'Filter29_thioester > 0'),
    ('CC(=S)S', 'Filter29_thioester > 0'),
    ('CC(=N)S', 'Filter29_thioester > 0'),
    ('CC(=O)CCBr', 'Filter26_alkyl_halide > 0'),
    ('CC(=O)CCI', 'Filter26_alkyl_halide > 0'),
    ('C(=O)CCI', 'Filter26_alkyl_halide > 0'),
    ('NC(=O)CCI', 'Filter26_alkyl_halide > 0'),
    ('CC=NOS(=O)N', 'Filter18_oxime_ester > 0'),
    ('NC=NOP(=O)N', 'Filter89_hydroxylamine > 0'),
    ('C=NOP(=O)N', 'Filter18_oxime_ester > 0'),
    ('SO', 'Filter31_so_bond > 0'),
    ('SOC', 'Filter31_so_bond > 0'),
    ('c1csocc1', 'OK'),
    ('OO', 'Filter32_oo_bond > 0'),
    ('COO', 'Filter32_oo_bond > 0'),
    ('c1coocc1', 'OK'),
    ('CC(=O)C=C', 'Filter44_michael_acceptor2 > 0'),
    ('C(=O)C=C', 'Filter38_aldehyde > 0'),
    ('OC(=O)C=C', 'Filter44_michael_acceptor2 > 0'),
    ('C1(=O)C=CC(=O)C=C1', 'Filter53_para_quinones > 0'),
    ('CIC', 'Filter49_halogen > 0'),
    ('CI(C)C', 'Filter49_halogen > 0'),
    ('CC=NOS(=O)N', 'Filter18_oxime_ester > 0'),
    ('NC=NO', 'Filter89_hydroxylamine > 0'),
    ('C(=O)NOS', 'Filter31_so_bond > 0'),
]
         
def test_hydrogen_suppression():
    alert_file_name = pkg_resources.resource_filename('rd_filters',
                                                      "data/alert_collection.csv")
    for rule_list, tests in [(["Inpharmatica"], test_inpharmatica),
                             (["LINT"], test_lint)]:
        rf = rd_filters.RDFilters(alert_file_name)
        rf.build_rule_list(rule_list)
        for smi, res in tests:
            result = rf.evaluate((smi,smi))
            assert result[2] == res, repr((result[:3], res))

