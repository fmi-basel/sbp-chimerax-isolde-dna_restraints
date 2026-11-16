"""Pure-Python (no ChimeraX imports) helpers for B-DNA restraints.

Everything here works on plain data and on "residue-like" objects that have
``.name``, ``.number``, ``.chain_id`` and ``.find_atom(name)``, so it can be
unit-tested outside ChimeraX (see ``tests/test_pure.py``).


* ``data/b_dna_dihedral_defs.json``    -- named DNA torsion definitions
  (was ``dictionaries/named_dihedrals.json`` -> ``b_dna``)
* ``data/b_dna_torsion_targets.json``  -- B-DNA target angles in degrees
  (was ``dictionaries/dna.json``)
* ``data/b_dna_distance_tables.json``  -- inter-/intra-nucleotide distance tables
  (were the ``nt_dist_inter_dict`` / ``nt_dist_intra_dict`` literals in
  ``restraints/restraint_utils.py: restrain_b_dna()``)
"""
import json
from collections import namedtuple
from functools import lru_cache
from math import radians
from pathlib import Path

_DATA_DIR = Path(__file__).resolve().parent / "data"

DNA_RESNAMES = ("DA", "DC", "DG", "DT")

# Watson-Crick partner of each nucleotide; the inter-chain distance tables are
# only valid for such pairs.
COMPLEMENT = {"DA": "DT", "DT": "DA", "DG": "DC", "DC": "DG"}

#: Torsion names restrained by ``isolde restrain b-dna``.
BACKBONE_TORSIONS = ("alpha", "beta", "gamma", "delta", "epsilon", "zeta")
SUGAR_TORSIONS = ("v0", "v1", "v2", "v3", "v4")
CHI_TORSIONS = ("chi1", "chi2")
TORSION_NAMES = BACKBONE_TORSIONS + SUGAR_TORSIONS + CHI_TORSIONS

DEFAULT_TORSION_SPRING_CONSTANT = 5000.0

#: Adaptive-distance-restraint parameters.  ``alpha`` is the same for all.
DISTANCE_ALPHA = -2
GROOVE_KAPPA = 100          # P-P groove-width restraints
INTER_KAPPA = 20            # base-pair (inter-chain) restraints
INTRA_KAPPA = 10            # within-nucleotide restraints
MINOR_GROOVE_P_P = 11.71    # Angstrom, +/- GROOVE_TOLERANCE
MAJOR_GROOVE_P_P = 17.94    # Angstrom, +/- GROOVE_TOLERANCE
GROOVE_TOLERANCE = 0.02
MINOR_GROOVE_OFFSET = 4     # residue i on chain 1 vs. i-4 on chain 2
MAJOR_GROOVE_OFFSET = 5     # residue i on chain 2 vs. i-5 on chain 1
MAX_WELL_HALF_WIDTH = 2.0   # c = min(target / 5, 2.0)


def well_half_width(target):
    """Adaptive-distance-restraint well half-width ``c`` for a given target."""
    return min(target / 5, MAX_WELL_HALF_WIDTH)


# ---------------------------------------------------------------------------
# data loading
# ---------------------------------------------------------------------------

@lru_cache(maxsize=None)
def _load(filename):
    with open(_DATA_DIR / filename, "rt") as f:
        return json.load(f)


def dihedral_definitions():
    """List of ``(resname, dihedral_name, atom_names, external_mask)`` ready to
    be passed to ``ProperDihedralMgr.add_dihedral_def()``.
    """
    d = _load("b_dna_dihedral_defs.json")
    out = []
    for resname in d["nucleotides"]:
        for dname, (atoms, mask) in d["bb_sugar"].items():
            out.append((resname, dname, list(atoms), list(mask)))
        for dname, (atoms, mask) in d["per_base"].get(resname, {}).items():
            out.append((resname, dname, list(atoms), list(mask)))
    return out


def torsion_targets():
    """``{resname: {torsion_name: target_in_radians}}`` for DA, DC, DG, DT.
    """
    t = _load("b_dna_torsion_targets.json")["b_dna_targets"]
    shared = {}
    for group in ("bb", "sugar"):
        for name, deg in t[group].items():
            shared[name] = radians(deg)
    out = {}
    for resname in DNA_RESNAMES:
        d = dict(shared)
        for name, deg in t[resname].items():
            d[name] = radians(deg)
        out[resname] = d
    return out


def distance_tables():
    """``{"inter": {resname: [[a1, a2, dist, tol], ...]}, "intra": {...}}``.

    For ``inter`` rows, ``a1`` is an atom of the residue on chain 1 and ``a2`` an
    atom of its Watson-Crick partner on chain 2.  For ``intra`` rows both atoms
    belong to the same residue.
    """
    return _load("b_dna_distance_tables.json")


# ---------------------------------------------------------------------------
# chain pairing
# ---------------------------------------------------------------------------

def pair_chains(residues):
    """Split ``residues`` into the two antiparallel DNA strands.

    Returns ``(strand1, strand2)``: ``strand1`` is sorted by ascending residue
    number and ``strand2`` by *descending* number, so that ``strand1[i]`` and
    ``strand2[i]`` form the i-th base pair.  The first strand is the one with
    the alphabetically first chain ID.

    Raises ``ValueError`` with a user-readable message on bad input.
    """
    by_chain = {}
    for r in residues:
        by_chain.setdefault(r.chain_id, []).append(r)
    if len(by_chain) != 2:
        raise ValueError("Expected residues from exactly two DNA chains, got %d (%s)"
                         % (len(by_chain), ", ".join(sorted(by_chain)) or "none"))
    c1, c2 = sorted(by_chain)
    strand1 = sorted(by_chain[c1], key=lambda r: r.number)
    strand2 = sorted(by_chain[c2], key=lambda r: r.number, reverse=True)
    if len(strand1) != len(strand2):
        raise ValueError("DNA chains must have the same number of residues "
                         "(chain %s: %d, chain %s: %d)"
                         % (c1, len(strand1), c2, len(strand2)))
    return strand1, strand2


def non_complementary_pairs(strand1, strand2):
    """Base pairs ``(res1, res2)`` that are not Watson-Crick (A-T / G-C)."""
    return [(a, b) for a, b in zip(strand1, strand2)
            if COMPLEMENT.get(a.name) != b.name]


# ---------------------------------------------------------------------------
# distance-restraint specifications
# ---------------------------------------------------------------------------

DistanceSpec = namedtuple(
    "DistanceSpec", "atom1 atom2 target tolerance kappa c alpha kind")


def build_distance_specs(strand1, strand2):
    """List of :class:`DistanceSpec` for ``isolde restrain b-dna``.

    ``strand1``/``strand2`` come from :func:`pair_chains`.  Pairs where either
    atom is missing (e.g. ``P`` of a 5'-terminal residue) are skipped.

    Order of the returned specs:
    minor-groove P-P, major-groove P-P, inter-chain base-pair distances, then
    intra-nucleotide distances.
    """
    tables = distance_tables()
    inter, intra = tables["inter"], tables["intra"]
    specs = []

    def add(atom1, atom2, target, tolerance, kappa, kind):
        if atom1 is None or atom2 is None:
            return
        specs.append(DistanceSpec(atom1, atom2, target, tolerance, kappa,
                                  well_half_width(target), DISTANCE_ALPHA, kind))

    # minor groove: residue i of strand 1 against residue i-4 of strand 2
    for i, res1 in enumerate(strand1):
        if i > 4:
            res2 = strand2[i - MINOR_GROOVE_OFFSET]
            add(res1.find_atom("P"), res2.find_atom("P"), MINOR_GROOVE_P_P,
                GROOVE_TOLERANCE, GROOVE_KAPPA, "minor_groove")

    # major groove: residue i of strand 2 against residue i-5 of strand 1
    for i, res2 in enumerate(strand2):
        if i > 4:
            res1 = strand1[i - MAJOR_GROOVE_OFFSET]
            add(res1.find_atom("P"), res2.find_atom("P"), MAJOR_GROOVE_P_P,
                GROOVE_TOLERANCE, GROOVE_KAPPA, "major_groove")

    # base pairs
    for res1, res2 in zip(strand1, strand2):
        for a1, a2, dist, tol in inter[res1.name]:
            add(res1.find_atom(a1), res2.find_atom(a2), dist, tol,
                INTER_KAPPA, "inter")

    # within each nucleotide
    for res1, res2 in zip(strand1, strand2):
        for res in (res1, res2):
            for a1, a2, dist, tol in intra[res.name]:
                add(res.find_atom(a1), res.find_atom(a2), dist, tol,
                    INTRA_KAPPA, "intra")
    return specs
