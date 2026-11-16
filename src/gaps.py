"""Pure-Python logic for ``isolde stepgap`` (no ChimeraX imports).

A *gap stop* is a residue at which the chain of residues "breaks":

* the first residue after a numbering gap (``number[j] - number[j-1] > 1``)
  within the same chain, and
* for a change of chain: the last residue of the old chain *and* the first
  residue of the new one.
"""
from bisect import bisect_left, bisect_right


def gap_stops(chain_ids, numbers):
    """Sorted list of indices that are gap stops.

    ``chain_ids`` and ``numbers`` are parallel sequences (one entry per residue,
    in stepping order).
    """
    if len(chain_ids) != len(numbers):
        raise ValueError("chain_ids and numbers must have the same length")
    stops = set()
    for j in range(1, len(chain_ids)):
        if chain_ids[j] != chain_ids[j - 1]:
            stops.add(j - 1)
            stops.add(j)
        elif numbers[j] - numbers[j - 1] > 1:
            stops.add(j)
    return sorted(stops)


def next_stop(stops, current):
    """First stop strictly after index ``current``, or ``None``."""
    k = bisect_right(stops, current)
    return stops[k] if k < len(stops) else None


def prev_stop(stops, current):
    """Last stop strictly before index ``current``, or ``None``."""
    k = bisect_left(stops, current)
    return stops[k - 1] if k > 0 else None
