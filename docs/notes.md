# Notes

RDKit was not able to parse some of the SMARTS in the structural_alerts table in ChEMBL.
I modified the SMARTS so that the patterns could be parsed.
I'm pretty sure I didn't alter the original intent of the SMARTS.
This document lists the changes I made and why I made them.

Each record consists of

- The rule number
- The original SMARTS
- The reason the RDKit couldn't parse the SMARTS
- My new version of the SMARTS

The file `alerts.csv` contains the edited SMARTS where applicable.

## rule 194

`[F,Cl,Br,I,$(O(S(=O)(=O)))]-[CH,CH2;!$(CF2)]-[N,n]`
As Andrew Dalke pointed out, `CF2` should actually be `C(F)F`
`[F,Cl,Br,I,$(O(S(=O)(=O)))]-[CH,CH2;!$(C(F)F)]-[N,n]`

## rule 196

`Parse Error: unclosed ring for input: '[N,n,O,S;!$(S(=O)(=O))]-[CH,CH2;!$(CF2)][F,Cl,Br,I,$(O(S(=O)(=O)))]'`
The RDKit isn't happy when `CF2` is not in brackets make it `[CF2]`
`[N,n,O,S;!$(S(=O)(=O))]-[CH,CH2;!$([CF2])][F,Cl,Br,I,$(O(S(=O)(=O)))]`

## rule 245

`SMARTS Parse Error: unclosed ring for input: '[CH2,$(CF2);R0][CH2,$(CF2);R0][CH2,$(CF2);R0][CH2,$(CF2);R0][CH2,$(CF2);R0][CH2,$(CF2);R0][CH2,$(CF2);R0][CH2,$(CF2);R0]'`
The RDKit isn't happy when `CF2` is not in brackets make it `[CF2]`
`[CH2,$([CF2]);R0][CH2,$([CF2]);R0][CH2,$([CF2]);R0][CH2,$([CF2]);R0][CH2,$([CF2]);R0][CH2,$([CF2]);R0][CH2,$([CF2]);R0][CH2,$([CF2]);R0]`

## rule 339

`[17:05:35] SMARTS Parse Error: syntax error for input: '[[CH;!R];!$(C-N)]=C([$(S(=O)(=O)),$(C(F)(F)(F)),$(C#N),$(N(=O)(=O)),$([N+](=O)[O-]),$(C(=O))])([$(S(=O)(=O)),$(C(F)(F)(F)),$(C#N),$(N(=O)(=O)),$([N+](=O)[O-]),$(C(=O))])'`
Remove unnecessary square brackets at the beginning on `[CH;!R]`
`[CH;!R;!$(C-N)]=C([$(S(=O)(=O)),$(C(F)(F)(F)),$(C#N),$(N(=O)(=O)),$([N+](=O)[O-]),$(C(=O))])([$(S(=O)(=O)),$(C(F)(F)(F)),$(C#N),$(N(=O)(=O)),$([N+](=O)[O-]),$(C(=O))])`

## rule 983

`SMARTS Parse Error: syntax error for input: '[2H,3H,11C,11c,14C,14c,125I,32P,33P,35S]'`
Versions of the RDKit prior to 2018-06 could not handle deuterium as `[2H]` or tritium as `[3H]` changed to `[2#1]` and `[3#1]`
`[2#1,3#1,11C,11c,14C,14c,125I,32P,33P,35S]`

## rule 1040

`SMARTS Parse Error: syntax error for input: '[N,C,S,O]-&!@[N,C,S,O]&!@[N,C,S,O]&!@[N,C,S,O]&!@[N,C,S,O]&!@[N,C,S,O]&!@[N,C,S,O]'`
The RDKit doesn't like implicit bond orders when specifying non-ring bonds e.g. `&!@` replace with `-&!@`
`[N,C,S,O]-&!@[N,C,S,O]-&!@[N,C,S,O]-&!@[N,C,S,O]-&!@[N,C,S,O]-&!@[N,C,S,O]-&!@[N,C,S,O]`

## rule 1053

`SMARTS Parse Error: syntax error for input: 'ac-*=&!@*-&!@C(=O)&!@ca'`
The RDKit doesn't like implicit bond orders when specifying non-ring bonds e.g. `&!@` replace with `-&!@`
`ac-*=&!@*-&!@C(=O)-&!@ca`

## rule 1072

`SMARTS Parse Error: syntax error for input: '[#6,#7]&!@[#6](=&!@[CH])&!@C(=O)-&!@[C,N,O,S]'`
The RDKit doesn't like implicit bond orders when specifying non-ring bonds e.g. `&!@` replace with `-&!@`
`[#6,#7]-&!@[#6](=&!@[CH])-&!@C(=O)-&!@[C,N,O,S]`

## rule 1079

`SMARTS Parse Error: syntax error for input: '*-C(=O)-&!@[NH]-C&!@C(=O)-&!@[NH]-*'`
The RDKit doesn't like implicit bond orders when specifying non-ring bonds e.g. `&!@` replace with `-&!@`
`*-C(=O)-&!@[NH]-C-&!@C(=O)-&!@[NH]-*`

## rule 1085

`SMARTS Parse Error: syntax error for input: '([#6]OP(=O)(*)O[#6].[#6]OP(=O)(*)O[#6].[#6]OP(=O)(*)O[#6])'`
Removed unnecessary leading and trailing parens
`[#6]OP(=O)(*)O[#6].[#6]OP(=O)(*)O[#6].[#6]OP(=O)(*)O[#6]`

## rule 1092

`SMARTS Parse Error: syntax error for input: 'c12cccc(C(=O)N(&!@C)C(=O)3)c2c3ccc1'`
The RDKit doesn't like implicit bond orders when specifying non-ring bonds e.g. `&!@` replace with `-&!@`
`c12cccc(C(=O)N(-&!@C)C(=O)3)c2c3ccc1`

## rule 1137

`[2H,3H,13C,14C,15N,125I,23F,22Na,32P,33P,35S,45Ca,57Co,103Ru,141Ce]`
Versions of the RDKit prior to 2018-06 could not handle deuterium as `[2H]` or tritium as `[3H]` changed to `[2#1]` and `[3#1]`
`[2#1,3#1,13C,14C,15N,125I,23F,22Na,32P,33P,35S,45Ca,57Co,103Ru,141Ce]`

## rule 1140

`SMARTS Parse Error: syntax error for input: 'Si~O'`
The RDKit can't handle silicon without brackets change to `[Si]`
`[Si]~O`
