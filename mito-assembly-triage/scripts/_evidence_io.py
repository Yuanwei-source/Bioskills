"""Input preservation and strict single-record DNA parsing for evidence tools."""
import os
from pathlib import Path
import re
import tempfile
import argparse
import sys


class InputArgumentParser(argparse.ArgumentParser):
    """Keep CLI input errors distinct from scientific REVIEW (exit 2)."""
    def error(self, message):
        self.print_usage(sys.stderr)
        self.exit(1, 'ERROR: %s\n' % message)


def read_single_fasta(path):
    name, pieces = None, []
    with open(path, encoding='utf-8') as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if name is not None:
                    raise ValueError('expected exactly one FASTA record: %s' % path)
                if not line[1:].strip():
                    raise ValueError('empty FASTA identifier: %s' % path)
                name = line[1:].split()[0]
            else:
                if name is None or not re.fullmatch('[ACGTRYSWKMBDHVNacgtryswkmbdhvn]+', line):
                    raise ValueError('invalid DNA FASTA sequence: %s' % path)
                pieces.append(line)
    if name is None or not pieces:
        raise ValueError('empty FASTA: %s' % path)
    return name, ''.join(pieces).upper()


def protect_output(output, inputs):
    output = Path(output)
    for source in inputs:
        source = Path(source)
        if output.resolve() == source.resolve() or (output.exists() and source.exists()
                                                     and os.path.samefile(output, source)):
            raise ValueError('output must not overwrite input: %s' % source)


def atomic_write_text(path, text):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temp = tempfile.mkstemp(prefix='.' + path.name + '.', dir=path.parent)
    try:
        with os.fdopen(fd, 'w', encoding='utf-8') as handle:
            handle.write(text)
        os.replace(temp, path)
    finally:
        if os.path.exists(temp):
            os.unlink(temp)
