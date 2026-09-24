#!/usr/bin/env python3
"""Shared GenBank/FASTA fixtures for the mito-assembly-triage regression tests.

Not a test module (no ``test_`` prefix) so ``unittest discover`` does not collect it.

Conventions
-----------
* ``write_gb`` feature specs use **0-based half-open** coordinates.
* ``segments`` (optional) builds a ``join()``/compound location; ``start``/``end``
  are then ignored.  This is how cross-origin genes are represented.
"""
import hashlib
import importlib.util
import random
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

#: Typical insect (Drosophila-type) arrangement, the anchor documented in
#: references/standard_gene_order.md.  Not a universal rule.
ORDER = [
    ("trnI", "tRNA", 70, "+"), ("trnQ", "tRNA", 68, "-"), ("trnM", "tRNA", 70, "+"),
    ("nad2", "CDS", 300, "+"), ("trnW", "tRNA", 69, "+"), ("trnC", "tRNA", 66, "-"),
    ("trnY", "tRNA", 68, "-"), ("cox1", "CDS", 300, "+"), ("trnL2", "tRNA", 70, "+"),
    ("cox2", "CDS", 300, "+"), ("trnK", "tRNA", 70, "+"), ("trnD", "tRNA", 68, "+"),
    ("atp8", "CDS", 300, "+"), ("atp6", "CDS", 300, "+"), ("cox3", "CDS", 300, "+"),
    ("trnG", "tRNA", 68, "+"), ("nad3", "CDS", 300, "+"), ("trnA", "tRNA", 70, "+"),
    ("trnR", "tRNA", 68, "+"), ("trnN", "tRNA", 70, "+"), ("trnS1", "tRNA", 68, "+"),
    ("trnE", "tRNA", 70, "+"), ("trnF", "tRNA", 67, "-"), ("nad5", "CDS", 300, "-"),
    ("trnH", "tRNA", 68, "-"), ("nad4", "CDS", 300, "-"), ("nad4l", "CDS", 300, "-"),
    ("trnT", "tRNA", 69, "+"), ("trnP", "tRNA", 68, "-"), ("nad6", "CDS", 300, "+"),
    ("cob", "CDS", 300, "+"), ("trnS2", "tRNA", 70, "+"), ("nad1", "CDS", 300, "-"),
    ("trnL1", "tRNA", 70, "-"), ("rrnL", "rRNA", 1300, "-"), ("trnV", "tRNA", 70, "-"),
    ("rrnS", "rRNA", 750, "-"),
]

#: table 5 sense codons usable inside a synthetic CDS (TAA/TAG are the stops;
#: AGA/AGG are Ser in table 5, so they are allowed here).
STOP_CODONS_TABLE5 = ("TAA", "TAG")
CODON_POOL = [
    a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT"
    if a + b + c not in STOP_CODONS_TABLE5
]


def load_module(name, relative_path):
    path = ROOT / relative_path
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def biopython_available():
    try:
        import Bio  # noqa: F401
    except ImportError:
        return False
    return True


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_gb(path, sequence, features, topology="circular", accession=None):
    """Write a GenBank record.

    features: list of dicts with keys
      name (str), type (str), strand ('+'/'-'),
      and either start/end (0-based half-open) or segments=[(s,e), ...],
      plus optional codon_start, anticodon, note, transl_except, product, gene.
    """
    from Bio.Seq import Seq
    from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
    from Bio.SeqRecord import SeqRecord
    from Bio import SeqIO

    seq_features = []
    for spec in features:
        if spec.get("segments"):
            parts = [FeatureLocation(a, b) for a, b in spec["segments"]]
            location = CompoundLocation(parts) if len(parts) > 1 else parts[0]
        else:
            location = FeatureLocation(spec["start"], spec["end"])
        strand = 1 if spec.get("strand", "+") == "+" else -1
        location = location._replace(strand=strand)
        qualifiers = {}
        if spec.get("gene"):
            qualifiers["gene"] = [spec["gene"]]
        if spec.get("product"):
            qualifiers["product"] = [spec["product"]]
        if spec.get("note"):
            qualifiers["note"] = [spec["note"]]
        if spec.get("codon_start"):
            qualifiers["codon_start"] = [str(spec["codon_start"])]
        if spec.get("anticodon"):
            qualifiers["anticodon"] = [spec["anticodon"]]
        if spec.get("transl_except"):
            qualifiers["transl_except"] = [spec["transl_except"]]
        seq_features.append(SeqFeature(location, type=spec["type"], qualifiers=qualifiers))

    annotations = {"molecule_type": "DNA"}
    if topology:
        annotations["topology"] = topology
    record = SeqRecord(Seq(sequence), id=accession or "test", name="test", description="",
                       annotations=annotations)
    record.features = seq_features
    SeqIO.write(record, path, "genbank")
    return Path(path)


def random_sequence(length, seed=7, skip_stops=False):
    rng = random.Random(seed)
    if skip_stops:
        return "".join(rng.choice(CODON_POOL) for _ in range(length // 3))
    return "".join(rng.choice("ACGT") for _ in range(length))


def cds_block(length=300, seed=1):
    """A stop-free CDS whose extracted (feature-orientation) sequence is ATG...TAA."""
    rng = random.Random(seed)
    return "ATG" + "".join(rng.choice(CODON_POOL) for _ in range((length - 6) // 3)) + "TAA"


def build_synthetic_ref(directory, order=ORDER, seed=7, split=3000, overlap=20):
    """Write ref.gb plus two overlapping scaffolds (s1.fa / s2.fa) under `directory`.

    The record is a positive control for the documented typical-insect arrangement:
    every CDS is stop-free in table 5 and on the expected strand.
    """
    from Bio.Seq import Seq

    blocks, features, position = [], [], 0
    for index, (name, kind, length, strand) in enumerate(order):
        if kind == "CDS":
            block = cds_block(length, seed=seed * 1000 + index)
            if strand == "-":
                block = str(Seq(block).reverse_complement())
        else:
            block = random_sequence(length, seed=seed * 1000 + index)
        assert len(block) == length, (name, len(block), length)
        blocks.append(block)
        features.append({"gene": name, "type": kind, "start": position,
                         "end": position + length, "strand": strand})
        position += length
    sequence = "".join(blocks)
    ref = write_gb(Path(directory) / "ref.gb", sequence, features)
    (Path(directory) / "s1.fa").write_text(">s1\n%s\n" % sequence[:split], encoding="utf-8")
    (Path(directory) / "s2.fa").write_text(">s2\n%s\n" % sequence[split - overlap:], encoding="utf-8")
    return ref, sequence


def run_script(relative_path, *args, timeout=600):
    return subprocess.run(
        [sys.executable, str(ROOT / relative_path), *[str(a) for a in args]],
        capture_output=True, text=True, timeout=timeout,
    )


def run_annot_check(path, *extra, timeout=120):
    return run_script(Path("scripts") / "annot_check.py", path, *extra, timeout=timeout)


def have(*tools):
    return all(shutil.which(tool) for tool in tools)
