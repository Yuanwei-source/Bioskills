"""Small SAM primitives. Coordinates returned here are 1-based, half-open.

SEQ/QUAL from SAM are already oriented to the reference (SAM v1 section 1.4).
Template identity is (read group, QNAME), not a claim of PCR independence.
"""
import re


def cigar_operations(cigar):
    if not re.fullmatch(r'(?:[1-9][0-9]*[MIDNSHP=X])+', cigar):
        raise ValueError('invalid mapped CIGAR: %r' % cigar)
    return [(int(n), op) for n, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar)]


def template_key(fields):
    if not fields[0] or fields[0] == '*':
        raise ValueError('QNAME is required for template counting')
    groups = [field[5:] for field in fields[11:] if field.startswith('RG:Z:')]
    if len(groups) > 1:
        raise ValueError('multiple RG tags')
    return (groups[0] if groups else '', fields[0])


def aligned_blocks(fields):
    """Contiguous M/= /X blocks; any insertion, deletion, skip or clip breaks one."""
    cursor = int(fields[3])
    if cursor < 1:
        raise ValueError('mapped POS must be positive')
    blocks, start = [], None
    for length, op in cigar_operations(fields[5]):
        if op in 'M=X':
            if start is None:
                start = cursor
            cursor += length
        else:
            if start is not None:
                blocks.append((start, cursor))
                start = None
            if op in 'DN':
                cursor += length
    if start is not None:
        blocks.append((start, cursor))
    return blocks
