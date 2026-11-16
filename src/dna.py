"""ChimeraX-side implementation of the B-DNA restraints.



Only public ISOLDE API is used (``ProperDihedralMgr.add_dihedral_def``,
``session_extensions.get_adaptive_dihedral_restraint_mgr``,
``restraint_utils.restrain_atom_pair_adaptive_distance`` ...), so ISOLDE itself
stays untouched.
"""
from . import dna_data

# Marker set on the session's ProperDihedralMgr once the DNA definitions have
# been added, so that they are added exactly once per manager instance.
_DEFS_FLAG = "_isolde_dna_defs_loaded"


def ensure_dna_dihedral_defs(session):
    """Register the named B-DNA dihedrals (alpha ... zeta, v0 ... v4, chi1, chi2)
    with ISOLDE's session-level ``ProperDihedralMgr``.  Safe to call repeatedly;
    returns the manager.

    ISOLDE's own manager only knows protein and glycan dihedrals.  The manager is
    a per-session singleton that is re-created for every new session, hence the
    flag is stored on the manager object itself.
    """
    from chimerax.isolde.session_extensions import get_proper_dihedral_mgr
    pdm = get_proper_dihedral_mgr(session)
    if getattr(pdm, _DEFS_FLAG, False):
        return pdm
    for resname, dname, atom_names, mask in dna_data.dihedral_definitions():
        pdm.add_dihedral_def(resname, dname, atom_names, mask)
    setattr(pdm, _DEFS_FLAG, True)
    return pdm


def _dna_residues(residues):
    """The subset of ``residues`` that are DA/DC/DG/DT."""
    import numpy
    return residues[numpy.isin(residues.names, dna_data.DNA_RESNAMES)]


def _single_structure(residues):
    from chimerax.core.errors import UserError
    us = residues.unique_structures
    if len(us) != 1:
        raise UserError("All residues must be from the same model!")
    return us[0]


# ---------------------------------------------------------------------------
# restrain
# ---------------------------------------------------------------------------

def restrain_b_dna_torsions(session, residues, spring_constant, reset=False):
    """Create adaptive dihedral restraints to B-DNA targets for every named DNA
    torsion in ``residues``.

    Existing restraints keep their targets unless ``reset`` is True, but all of them get the new ``spring_constant``.  Returns the
    total number of restraints now covering the selection.
    """
    from chimerax.isolde import session_extensions as sx
    ensure_dna_dihedral_defs(session)
    structure = _single_structure(residues)
    mgr = sx.get_adaptive_dihedral_restraint_mgr(structure)
    targets = dna_data.torsion_targets()
    resnames = residues.unique_names

    total = 0
    for torsion in dna_data.TORSION_NAMES:
        for rname in resnames:
            target = targets[rname][torsion]
            current = residues[residues.names == rname]
            existing = None
            if not reset:
                existing = mgr.get_restraints_by_residues_and_name(current, torsion)
            restraints = mgr.add_restraints_by_residues_and_name(current, torsion)
            if existing is not None:
                restraints = restraints.subtract(existing)
            if not len(restraints):
                continue
            restraints.targets = target
            restraints.displays = True
            restraints.enableds = True
        covered = mgr.get_restraints_by_residues_and_name(residues, torsion)
        covered.spring_constants = spring_constant
        total += len(covered)
    return total


def restrain_b_dna_distances(residues):
    """Adaptive distance restraints that hold two paired DNA strands in B-form
    (groove-width P-P restraints, base-pair distances, intra-nucleotide
    distances).  Returns ``(n_created, n_skipped)``.
    """
    from chimerax.core.errors import UserError
    from chimerax.isolde.restraints.restraint_utils import \
        restrain_atom_pair_adaptive_distance
    try:
        strand1, strand2 = dna_data.pair_chains(list(residues))
    except ValueError as e:
        raise UserError(str(e))

    log = residues.unique_structures[0].session.logger
    mismatched = dna_data.non_complementary_pairs(strand1, strand2)
    if mismatched:
        log.warning(
            "isolde restrain b-dna: %d base pair(s) are not Watson-Crick "
            "complementary (first: %s - %s). The base-pair distance tables "
            "assume A-T/G-C pairing, so the pairing of the two chains may be "
            "wrong (chains are paired first-to-last)."
            % (len(mismatched), mismatched[0][0], mismatched[0][1]))

    specs = dna_data.build_distance_specs(strand1, strand2)
    created = skipped = 0
    for s in specs:
        try:
            restrain_atom_pair_adaptive_distance(
                s.atom1, s.atom2, s.target, s.tolerance, s.kappa, s.c, s.alpha)
            created += 1
        except ValueError:
            # e.g. atoms that cannot be restrained against each other
            skipped += 1
    return created, skipped


def restrain_b_dna(session, residues, spring_constant, torsions=True,
                   distances=True, reset=False):
    """Entry point for ``isolde restrain b-dna``."""
    from chimerax.core.errors import UserError
    if not len(residues):
        raise UserError("isolde restrain b-dna: no residues specified.")
    bad = sorted(set(residues.names) - set(dna_data.DNA_RESNAMES))
    if bad:
        raise UserError(
            "isolde restrain b-dna: only DNA residues (DA, DC, DG, DT) can be "
            "restrained; the selection also contains: %s" % ", ".join(bad))
    _single_structure(residues)
    log = session.logger
    if torsions:
        n = restrain_b_dna_torsions(session, residues, spring_constant, reset=reset)
        log.info("isolde restrain b-dna: %d torsion restraints on %d residues "
                 "(spring constant %g kJ/mol/rad^2)" % (n, len(residues), spring_constant))
    if distances:
        created, skipped = restrain_b_dna_distances(residues)
        msg = "isolde restrain b-dna: %d adaptive distance restraints created" % created
        if skipped:
            msg += " (%d pairs skipped)" % skipped
        log.info(msg)


# ---------------------------------------------------------------------------
# adjust / release
# ---------------------------------------------------------------------------

def _torsion_names(backbone, sidechains):
    names = []
    if backbone:
        names.extend(dna_data.BACKBONE_TORSIONS + dna_data.SUGAR_TORSIONS)
    if sidechains:
        names.extend(dna_data.CHI_TORSIONS)
    return names


def adjust_b_dna_torsions(session, residues, backbone=True, sidechains=True,
                          angle_range=None, spring_constant=None, alpha=None):
    """Change kappa (via ``angle_range``, degrees), spring constant and/or alpha
    of existing B-DNA torsion restraints.  Returns the number of restraints
    touched."""
    from math import radians
    residues = _dna_residues(residues)
    if not len(residues):
        return 0
    ensure_dna_dihedral_defs(session)
    from chimerax.isolde import session_extensions as sx
    kappa = None
    if angle_range is not None:
        from chimerax.isolde.molobject import AdaptiveDihedralRestraintMgr
        kappa = AdaptiveDihedralRestraintMgr.angle_range_to_kappa(radians(angle_range))
    n = 0
    for m in residues.unique_structures:
        mgr = sx.get_adaptive_dihedral_restraint_mgr(m, create=False)
        if mgr is None:
            continue
        m_res = residues[residues.structures == m]
        for name in _torsion_names(backbone, sidechains):
            restraints = mgr.get_restraints_by_residues_and_name(m_res, name)
            if not len(restraints):
                continue
            if kappa is not None:
                restraints.kappas = kappa
            if spring_constant is not None:
                restraints.spring_constants = spring_constant
            if alpha is not None:
                restraints.alphas = alpha
            n += len(restraints)
    return n


def release_b_dna_torsions(session, residues, backbone=True, sidechains=True):
    """Disable existing B-DNA torsion restraints.  Returns the number released."""
    residues = _dna_residues(residues)
    if not len(residues):
        return 0
    ensure_dna_dihedral_defs(session)
    from chimerax.isolde import session_extensions as sx
    n = 0
    for m in residues.unique_structures:
        mgr = sx.get_adaptive_dihedral_restraint_mgr(m, create=False)
        if mgr is None:
            continue
        m_res = residues[residues.structures == m]
        for name in _torsion_names(backbone, sidechains):
            restraints = mgr.get_restraints_by_residues_and_name(m_res, name)
            restraints.enableds = False
            n += len(restraints)
    return n
