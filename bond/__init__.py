"""BOND: Bond-centric Omics Nexus Decipherer.

Public API for bond-type annotation of metabolite structures.

>>> from bond import BondJudge
>>> judge = BondJudge()
>>> result = judge.judge("OC[C@H]1O[C@@H](Oc2ccccc2)[C@H](O)[C@@H](O)[C@@H]1O")
>>> sorted(result.bond_codes)
['G14']
"""

from bond.judge_core import (  # noqa: F401
    BondJudge,
    JudgeResult,
    detect_glycoside,
    is_sugar_ring,
    run_table,
    scan_nongly,
)

__version__ = "0.1.0"
