"""Commands added to ISOLDE by the ChimeraX-ISOLDE-DNA_RESTRAINTS bundle.

``isolde restrain b-dna``   restrain DNA torsions + inter/intra distances to B-form
``isolde adjust b-dna``     change strength / width of those torsion restraints
``isolde release b-dna``    switch those torsion restraints off
``isolde stepgap``          jump to the next/previous numbering gap or chain break
"""
from chimerax.core.commands import (
    CmdDesc, register, BoolArg, FloatArg, PositiveIntArg, EnumOf,
)
from chimerax.core.errors import UserError

_HELP_URL = "help:user/commands/isolde_dna_restraints.html"


def register_command(command_name, logger):
    """BundleAPI dispatch: one registration per declared command name."""
    from chimerax.atomic import ResiduesArg, AtomicStructureArg
    if command_name == "isolde restrain b-dna":
        desc = CmdDesc(
            required=[("residues", ResiduesArg)],
            keyword=[
                ("spring_constant", FloatArg),
                ("torsions", BoolArg),
                ("distances", BoolArg),
                ("reset", BoolArg),
            ],
            synopsis="Restrain two paired DNA chains to B-form (EXPERIMENTAL)",
            url=_HELP_URL + "#restrain-b-dna",
        )
        register(command_name, desc, restrain_b_dna_cmd, logger=logger)
    elif command_name == "isolde adjust b-dna":
        desc = CmdDesc(
            required=[("residues", ResiduesArg)],
            keyword=[
                ("backbone", BoolArg),
                ("sidechains", BoolArg),
                ("angle_range", FloatArg),
                ("alpha", FloatArg),
                ("spring_constant", FloatArg),
            ],
            synopsis="Adjust existing B-DNA torsion restraints",
            url=_HELP_URL + "#adjust-b-dna",
        )
        register(command_name, desc, adjust_b_dna_cmd, logger=logger)
    elif command_name == "isolde release b-dna":
        desc = CmdDesc(
            required=[("residues", ResiduesArg)],
            keyword=[
                ("backbone", BoolArg),
                ("sidechains", BoolArg),
            ],
            synopsis="Release B-DNA torsion restraints",
            url=_HELP_URL + "#release-b-dna",
        )
        register(command_name, desc, release_b_dna_cmd, logger=logger)
    elif command_name == "isolde stepgap":
        desc = CmdDesc(
            optional=[("direction", EnumOf(("next", "prev")))],
            keyword=[
                ("model", AtomicStructureArg),
                ("polymeric_only", BoolArg),
                ("select", BoolArg),
                ("view_distance", FloatArg),
                ("interpolate_frames", PositiveIntArg),
            ],
            synopsis="Step to the next/previous numbering gap or chain break",
            url=_HELP_URL + "#stepgap",
        )
        register(command_name, desc, stepgap_cmd, logger=logger)


# ---------------------------------------------------------------------------
# isolde restrain / adjust / release b-dna
# ---------------------------------------------------------------------------

def restrain_b_dna_cmd(session, residues, spring_constant=None, torsions=True,
                       distances=True, reset=False):
    from . import dna, dna_data
    if spring_constant is None:
        spring_constant = dna_data.DEFAULT_TORSION_SPRING_CONSTANT
    dna.restrain_b_dna(session, residues, spring_constant, torsions=torsions,
                       distances=distances, reset=reset)


def adjust_b_dna_cmd(session, residues, backbone=True, sidechains=True,
                     angle_range=None, alpha=None, spring_constant=None):
    from . import dna
    n = dna.adjust_b_dna_torsions(
        session, residues, backbone=backbone, sidechains=sidechains,
        angle_range=angle_range, spring_constant=spring_constant, alpha=alpha)
    session.logger.info("isolde adjust b-dna: adjusted %d torsion restraints" % n)


def release_b_dna_cmd(session, residues, backbone=True, sidechains=True):
    from . import dna
    n = dna.release_b_dna_torsions(
        session, residues, backbone=backbone, sidechains=sidechains)
    session.logger.info("isolde release b-dna: released %d torsion restraints" % n)


# ---------------------------------------------------------------------------
# isolde stepgap
# ---------------------------------------------------------------------------

def _pick_model(session, model):
    if model is not None:
        return model
    m = None
    isolde = getattr(session, "isolde", None)
    if isolde is not None:
        m = isolde.selected_model
    if m is None:
        from chimerax.atomic import AtomicStructure
        structures = session.models.list(type=AtomicStructure)
        if structures:
            m = structures[0]
    if m is None:
        raise UserError("No atomic structures are open!")
    return m


def stepgap_cmd(session, direction="next", model=None, polymeric_only=True,
                select=True, view_distance=None, interpolate_frames=None):
    """Move ISOLDE's residue stepper to the next/previous gap stop.

    See :mod:`.gaps` for what counts as a gap stop.
    """
    from chimerax.atomic import Residue
    from chimerax.isolde.navigate import get_stepper
    from . import gaps

    m = _pick_model(session, model)
    residues = m.residues
    if polymeric_only:
        polymeric = residues[residues.polymer_types != Residue.PT_NONE]
        if len(polymeric):
            residues = polymeric
        else:
            session.logger.warning(
                "isolde stepgap: this model contains no polymeric residues; "
                "ignoring polymericOnly.")
    if len(residues) < 2:
        raise UserError("isolde stepgap: need at least two residues.")

    stops = gaps.gap_stops(list(residues.chain_ids), [int(n) for n in residues.numbers])
    if not stops:
        session.logger.info("isolde stepgap: no numbering gaps or chain breaks in %s" % m)
        return

    rs = get_stepper(m)
    if view_distance is not None:
        rs.view_distance = view_distance
    if interpolate_frames is not None:
        rs.interpolate_frames = interpolate_frames

    cr = rs.current_residue
    ci = residues.index(cr) if cr is not None else -1
    if ci == -1:
        # no usable current residue: next starts from the first residue,
        # prev from the end of the model
        ci = 0 if direction == "next" else len(residues)
    target = (gaps.next_stop if direction == "next" else gaps.prev_stop)(stops, ci)
    if target is None:
        session.logger.info("isolde stepgap: no further gaps in the '%s' direction"
                            % direction)
        return
    if select:
        session.selection.clear()
    res = residues[target]
    rs.step_to(res)
    if select:
        atoms = res.atoms
        atoms.selected = True
        atoms.intra_bonds.selected = True
    session.logger.info("isolde stepgap: stepped to %s" % res)
