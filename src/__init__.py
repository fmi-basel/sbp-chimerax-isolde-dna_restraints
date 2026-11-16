# ChimeraX-ISOLDE-DNA_RESTRAINTS: Extensions to ISOLDE
# Installed as the package chimerax.isolde.dna_restraints (see pyproject.toml).

__version__ = "0.1.0"

from chimerax.core.toolshed import BundleAPI


class _DNARESTRAINTSAPI(BundleAPI):

    api_version = 1

    @staticmethod
    def register_command(bi, ci, logger):
        # Called lazily, once per command name declared under
        # [tool.chimerax.command."..."] in pyproject.toml, the first time
        # that command is used.
        from . import cmd
        cmd.register_command(ci.name, logger)


bundle_api = _DNARESTRAINTSAPI()
