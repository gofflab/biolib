"""Sequence data structures and molecular biology utilities.

Provides SeqDict, a dictionary subclass for ordered molecular sequences, and
a variety of constants and functions for DNA/RNA/protein operations including
codon translation, reverse complementation, GC content calculation, and
Kimura sequence evolution simulation.

Author: lgoff (derived from rasmus seqlib)
"""

import copy
import math
import random

# from rasmus import util  # NOTE: rasmus is not available; util functions inlined below


class SeqDict (dict):
    """A dictionary for molecular sequences that also tracks insertion order.

    Useful for reading and writing sequences from FASTA files where order
    matters. Keys are sequence names; values are sequence strings. See
    fasta.FastaDict for a subclass that implements FASTA reading and writing.

    Attributes:
        names: List of sequence names in insertion order.
    """

    def __init__(self):
        """Initialize an empty SeqDict."""
        dict.__init__(self)

        self.names = []


    def orderNames(self, aln):
        """Reorder self.names to match the key order of another dict.

        Args:
            aln: A dict (typically another SeqDict or alignment) whose key
                order is used to sort self.names.
        """

        # Inlined util.list2lookup: creates a dict mapping list items to their index
        lookup = {v: i for i, v in enumerate(aln.keys())}
        self.names.sort(key=lambda x: lookup[x])


    # add a key, value pair
    def add(self, key, value, errors=False):
        """Add a key-value pair, keeping the longest value on duplicate keys.

        If the key already exists and the new value is at least as long as the
        stored value, the stored value is replaced. The insertion order in
        self.names is preserved (duplicate keys do not add to names).

        Args:
            key: Sequence name string.
            value: Sequence string.
            errors: If True, write a warning to stderr on duplicate keys
                (default False).
        """
        if key in self:
            if errors:
                # Inlined util.logger: write to stderr
                import sys
                sys.stderr.write("duplicate key %s\n" % str(key))

            # keep the longest value, by default
            if len(value) >= len(self[key]):
                dict.__setitem__(self, key, value)
        else:
            self.names.append(key)
            dict.__setitem__(self, key, value)


    def get(self, keys, new=None):
        """Return a new SeqDict containing only the given keys.

        Args:
            keys: Iterable of key names to include.
            new: Optional pre-existing SeqDict to populate. If None, a new
                instance of the same type is created.

        Returns:
            A SeqDict (or instance of the same subclass) containing the
            requested keys that are present in self.
        """

        if new == None:
            new = type(self)()

        for key in keys:
            if key in self:
                new[key] = self[key]

        return new


    def alignlen(self):
        """
        If this SeqDict is an alignment, this function
        will return its length
        """

        return len(list(self.values())[0])


    # The following methods keep names in sync with dictionary keys
    def __setitem__(self, key, value):
        """Set a key-value pair and add key to self.names if new."""
        if key not in self:
            self.names.append(key)
        dict.__setitem__(self, key, value)

    def __delitem__(self, key):
        """Delete a key and remove it from self.names."""
        self.names.remove(key)

    def update(self, dct):
        """Update from another dict, appending new keys to self.names.

        Args:
            dct: Dict-like object whose items will be merged into self.
        """
        for key in dct:
            if key not in self.names:
                self.names.append(key)
        dict.update(self, dct)

    def setdefault(self, key, value):
        """Set key to value only if key is absent, tracking order.

        Args:
            key: Key to look up or set.
            value: Default value to assign if key is missing.
        """
        if key not in self.names:
            self.names.append(key)
        dict.setdefault(self, key, value)

    def clear(self):
        """Remove all items and reset self.names to an empty list."""
        self.names = []
        dict.clear(self)

    # keys are always sorted in order added
    def keys(self):
        """Return keys in insertion order.

        Returns:
            List of key names in insertion order.
        """
        return list(self.names)

    def iterkeys(self):
        """Iterate over keys in insertion order.

        Returns:
            Iterator over key name strings.
        """
        return iter(self.names)

    def values(self):
        """Return values in key insertion order.

        Returns:
            List of sequence strings in the same order as self.names.
        """
        return [self[key] for key in self.iterkeys()]

    def itervalues(self):
        """Iterate over values in key insertion order.

        Returns:
            Generator yielding sequence strings in insertion order.
        """
        def func():
            for key in self.iterkeys():
                yield self[key]
        return func()

    def iteritems(self):
        """Iterate over (key, value) pairs in key insertion order.

        Returns:
            Generator yielding (name, sequence) tuples.
        """
        def func():
            for key in self.iterkeys():
                yield (key, self[key])
        return func()

    def items(self):
        """Return list of (key, value) pairs in insertion order.

        Returns:
            List of (name, sequence) tuples.
        """
        return list(self.iteritems())

    def __iter__(self):
        """Iterate over keys in insertion order."""
        return iter(self.names)

    def __len__(self):
        """Return the number of sequences stored."""
        return len(self.names)



#--------------------------------------------------------------------------------
# Constants
#--------------------------------------------------------------------------------


# standard codon table
CODON_TABLE = {
    "TTT": "F",  "CTT": "L",  "ATT": "I",  "GTT": "V",
    "TTC": "F",  "CTC": "L",  "ATC": "I",  "GTC": "V",
    "TTA": "L",  "CTA": "L",  "ATA": "I",  "GTA": "V",
    "TTG": "L",  "CTG": "L",  "ATG": "M",  "GTG": "V",

    "TCT": "S",  "CCT": "P",  "ACT": "T",  "GCT": "A",
    "TCC": "S",  "CCC": "P",  "ACC": "T",  "GCC": "A",
    "TCA": "S",  "CCA": "P",  "ACA": "T",  "GCA": "A",
    "TCG": "S",  "CCG": "P",  "ACG": "T",  "GCG": "A",

    "TAT": "Y",  "CAT": "H",  "AAT": "N",  "GAT": "D",
    "TAC": "Y",  "CAC": "H",  "AAC": "N",  "GAC": "D",
    "TAA": "*",  "CAA": "Q",  "AAA": "K",  "GAA": "E",
    "TAG": "*",  "CAG": "Q",  "AAG": "K",  "GAG": "E",

    "TGT": "C",  "CGT": "R",  "AGT": "S",  "GGT": "G",
    "TGC": "C",  "CGC": "R",  "AGC": "S",  "GGC": "G",
    "TGA": "*",  "CGA": "R",  "AGA": "R",  "GGA": "G",
    "TGG": "W",  "CGG": "R",  "AGG": "R",  "GGG": "G",

    "---": "-"
}

# codon table specific to the Candida species
CANDIDA_CODON_TABLE = copy.copy(CODON_TABLE)
CANDIDA_CODON_TABLE["CTG"] = "S"  # originally L


# make reverse codon table
REV_CODON_TABLE = {}
for key,val in CODON_TABLE.items():
    REV_CODON_TABLE.setdefault(val, []).append(key)


# make degenerate counts
#
# example:
#
# CGT => "R"
# CGC => "R"
# CGA => "R"
# CGG => "R"
#
# CODON_DEGEN["R"] = [1, 1, 4]
# CODON_DEGEN["CGT"] = [1, 1, 4]
#
CODON_DEGEN = {}
AA_DEGEN = {}
for aa, lst in REV_CODON_TABLE.items():
    # Inlined: map(lambda x: len(util.unique(x)), zip(*lst))
    # util.unique(x) returns unique elements; replaced with set(x)
    folds = [len(set(x)) for x in zip(* lst)]
    for codon in lst:
        AA_DEGEN[aa] = folds
        CODON_DEGEN[codon] = folds


# substitution types
SUB_NONE = 0  # none
SUB_TSIT = 1  # tranSition
SUB_TVER = 2  # transVersion
SUB_INS  = 3  # insert
SUB_DEL  = 4  # del
SUBSTITUTION_TYPES = {
    "AA": SUB_NONE, "AC": SUB_TVER, "AG": SUB_TSIT, "AT": SUB_TVER,
    "CA": SUB_TVER, "CC": SUB_NONE, "CG": SUB_TVER, "CT": SUB_TSIT,
    "GA": SUB_TSIT, "GC": SUB_TVER, "GG": SUB_NONE, "GT": SUB_TVER,
    "TA": SUB_TVER, "TC": SUB_TSIT, "TG": SUB_TVER, "TT": SUB_NONE,

    "A-": SUB_DEL, "C-": SUB_DEL, "G-": SUB_DEL, "T-": SUB_DEL,
    "-A": SUB_INS, "-C": SUB_INS, "-G": SUB_INS, "-T": SUB_INS,

    "--": SUB_NONE, "NN": SUB_NONE,
    "NA": SUB_NONE, "NC": SUB_NONE, "NT": SUB_NONE, "NG": SUB_NONE,
    "AN": SUB_NONE, "CN": SUB_NONE, "TN": SUB_NONE, "GN": SUB_NONE,
    "N-": SUB_NONE, "N-": SUB_NONE, "N-": SUB_NONE, "N-": SUB_NONE,
    "-N": SUB_NONE, "-N": SUB_NONE, "-N": SUB_NONE, "-N": SUB_NONE
}


# hydrophobic / hydrophilic
def hydrophobic(aa):
    """Return a numeric hydrophobicity score for a single amino acid.

    Args:
        aa: Single-letter amino-acid code string.

    Returns:
        2.0 for strongly hydrophobic residues (VILMFWC),
        1.0 for weakly hydrophobic residues (AYHTSPG),
        0.5 for weakly hydrophilic residues (RK),
        0.0 for all other residues.
    """
    if aa in 'VILMFWC': return 2.0
    if aa in 'AYHTSPG': return 1.0
    if aa in 'RK': return 0.5
    return 0.0


AA_PROPERTY = {'A': 'weakly hydrophobic',
               'R': 'charged',
               'N': 'polar',
               'D': 'charged',
               'C': 'polar',
               'E': 'charged',
               'Q': 'polar',
               'G': 'turn',
               'H': 'charged',
               'I': 'hydrophobic',
               'L': 'hydrophobic',
               'K': 'polar',
               'M': 'met',
               'F': 'hydrophobic',
               'P': 'hydrophobic',
               'S': 'polar',
               'T': 'polar',
               'W': 'hydrophobic',
               'Y': 'polar',
               'V': 'hydrophobic',
               'U': 'polar',
               '*': 'stop',
               '-': 'gap'}



BLOSUM62 = \
       {'A': {'A': 4, 'R':-1, 'N':-2, 'D':-2, 'C': 0, 'Q':-1, 'E':-1, 'G': 0, 'H':-2, 'I':-1, 'L':-1, 'K':-1,
              'M':-1, 'F':-2, 'P':-1, 'S': 1, 'T': 0, 'W':-3, 'Y':-2, 'V': 0, 'B':-2, 'Z':-1, 'X': 0, '*':-4},
        'R': {'A':-1, 'R': 5, 'N': 0, 'D':-2, 'C':-3, 'Q': 1, 'E': 0, 'G':-2, 'H': 0, 'I':-3, 'L':-2, 'K': 2,
              'M':-1, 'F':-3, 'P':-2, 'S':-1, 'T':-1, 'W':-3, 'Y':-2, 'V':-3, 'B':-1, 'Z': 0, 'X':-1, '*':-4},
        'N': {'A':-2, 'R': 0, 'N': 6, 'D': 1, 'C':-3, 'Q': 0, 'E': 0, 'G': 0, 'H': 1, 'I':-3, 'L':-3, 'K': 0,
              'M':-2, 'F':-3, 'P':-2, 'S': 1, 'T': 0, 'W':-4, 'Y':-2, 'V':-3, 'B': 3, 'Z': 0, 'X':-1, '*':-4},
        'D': {'A':-2, 'R':-2, 'N': 1, 'D': 6, 'C':-3, 'Q': 0, 'E': 2, 'G':-1, 'H':-1, 'I':-3, 'L':-4, 'K':-1,
              'M':-3, 'F':-3, 'P':-1, 'S': 0, 'T':-1, 'W':-4, 'Y':-3, 'V':-3, 'B': 4, 'Z': 1, 'X':-1, '*':-4},
        'C': {'A': 0, 'R':-3, 'N':-3, 'D':-3, 'C': 9, 'Q':-3, 'E':-4, 'G':-3, 'H':-3, 'I':-1, 'L':-1, 'K':-3,
              'M':-1, 'F':-2, 'P':-3, 'S':-1, 'T':-1, 'W':-2, 'Y':-2, 'V':-1, 'B':-3, 'Z':-3, 'X':-2, '*':-4},
        'Q': {'A':-1, 'R': 1, 'N': 0, 'D': 0, 'C':-3, 'Q': 5, 'E': 2, 'G':-2, 'H': 0, 'I':-3, 'L':-2, 'K': 1,
              'M': 0, 'F':-3, 'P':-1, 'S': 0, 'T':-1, 'W':-2, 'Y':-1, 'V':-2, 'B': 0, 'Z': 3, 'X':-1, '*':-4},
        'E': {'A':-1, 'R': 0, 'N': 0, 'D': 2, 'C':-4, 'Q': 2, 'E': 5, 'G':-2, 'H': 0, 'I':-3, 'L':-3, 'K': 1,
              'M':-2, 'F':-3, 'P':-1, 'S': 0, 'T':-1, 'W':-3, 'Y':-2, 'V':-2, 'B': 1, 'Z': 4, 'X':-1, '*':-4},
        'G': {'A': 0, 'R':-2, 'N': 0, 'D':-1, 'C':-3, 'Q':-2, 'E':-2, 'G': 6, 'H':-2, 'I':-4, 'L':-4, 'K':-2,
              'M':-3, 'F':-3, 'P':-2, 'S': 0, 'T':-2, 'W':-2, 'Y':-3, 'V':-3, 'B':-1, 'Z':-2, 'X':-1, '*':-4},
        'H': {'A':-2, 'R': 0, 'N': 1, 'D':-1, 'C':-3, 'Q': 0, 'E': 0, 'G':-2, 'H': 8, 'I':-3, 'L':-3, 'K':-1,
              'M':-2, 'F':-1, 'P':-2, 'S':-1, 'T':-2, 'W':-2, 'Y': 2, 'V':-3, 'B': 0, 'Z': 0, 'X':-1, '*':-4},
        'I': {'A':-1, 'R':-3, 'N':-3, 'D':-3, 'C':-1, 'Q':-3, 'E':-3, 'G':-4, 'H':-3, 'I': 4, 'L': 2, 'K':-3,
              'M': 1, 'F': 0, 'P':-3, 'S':-2, 'T':-1, 'W':-3, 'Y':-1, 'V': 3, 'B':-3, 'Z':-3, 'X':-1, '*':-4},
        'L': {'A':-1, 'R':-2, 'N':-3, 'D':-4, 'C':-1, 'Q':-2, 'E':-3, 'G':-4, 'H':-3, 'I': 2, 'L': 4, 'K':-2,
              'M': 2, 'F': 0, 'P':-3, 'S':-2, 'T':-1, 'W':-2, 'Y':-1, 'V': 1, 'B':-4, 'Z':-3, 'X':-1, '*':-4},
        'K': {'A':-1, 'R': 2, 'N': 0, 'D':-1, 'C':-3, 'Q': 1, 'E': 1, 'G':-2, 'H':-1, 'I':-3, 'L':-2, 'K': 5,
              'M':-1, 'F':-3, 'P':-1, 'S': 0, 'T':-1, 'W':-3, 'Y':-2, 'V':-2, 'B': 0, 'Z': 1, 'X':-1, '*':-4},
        'M': {'A':-1, 'R':-1, 'N':-2, 'D':-3, 'C':-1, 'Q': 0, 'E':-2, 'G':-3, 'H':-2, 'I': 1, 'L': 2, 'K':-1,
              'M': 5, 'F': 0, 'P':-2, 'S':-1, 'T':-1, 'W':-1, 'Y':-1, 'V': 1, 'B':-3, 'Z':-1, 'X':-1, '*':-4},
        'F': {'A':-2, 'R':-3, 'N':-3, 'D':-3, 'C':-2, 'Q':-3, 'E':-3, 'G':-3, 'H':-1, 'I': 0, 'L': 0, 'K':-3,
              'M': 0, 'F': 6, 'P':-4, 'S':-2, 'T':-2, 'W': 1, 'Y': 3, 'V':-1, 'B':-3, 'Z':-3, 'X':-1, '*':-4},
        'P': {'A':-1, 'R':-2, 'N':-2, 'D':-1, 'C':-3, 'Q':-1, 'E':-1, 'G':-2, 'H':-2, 'I':-3, 'L':-3, 'K':-1,
              'M':-2, 'F':-4, 'P': 7, 'S':-1, 'T':-1, 'W':-4, 'Y':-3, 'V':-2, 'B':-2, 'Z':-1, 'X':-2, '*':-4},
        'S': {'A': 1, 'R':-1, 'N': 1, 'D': 0, 'C':-1, 'Q': 0, 'E': 0, 'G': 0, 'H':-1, 'I':-2, 'L':-2, 'K': 0,
              'M':-1, 'F':-2, 'P':-1, 'S': 4, 'T': 1, 'W':-3, 'Y':-2, 'V':-2, 'B': 0, 'Z': 0, 'X': 0, '*':-4},
        'T': {'A': 0, 'R':-1, 'N': 0, 'D':-1, 'C':-1, 'Q':-1, 'E':-1, 'G':-2, 'H':-2, 'I':-1, 'L':-1, 'K':-1,
              'M':-1, 'F':-2, 'P':-1, 'S': 1, 'T': 5, 'W':-2, 'Y':-2, 'V': 0, 'B':-1, 'Z':-1, 'X': 0, '*':-4},
        'W': {'A':-3, 'R':-3, 'N':-4, 'D':-4, 'C':-2, 'Q':-2, 'E':-3, 'G':-2, 'H':-2, 'I':-3, 'L':-2, 'K':-3,
              'M':-1, 'F': 1, 'P':-4, 'S':-3, 'T':-2, 'W':11, 'Y': 2, 'V':-3, 'B':-4, 'Z':-3, 'X':-2, '*':-4},
        'Y': {'A':-2, 'R':-2, 'N':-2, 'D':-3, 'C':-2, 'Q':-1, 'E':-2, 'G':-3, 'H': 2, 'I':-1, 'L':-1, 'K':-2,
              'M':-1, 'F': 3, 'P':-3, 'S':-2, 'T':-2, 'W': 2, 'Y': 7, 'V':-1, 'B':-3, 'Z':-2, 'X':-1, '*':-4},
        'V': {'A': 0, 'R':-3, 'N':-3, 'D':-3, 'C':-1, 'Q':-2, 'E':-2, 'G':-3, 'H':-3, 'I': 3, 'L': 1, 'K':-2,
              'M': 1, 'F':-1, 'P':-2, 'S':-2, 'T': 0, 'W':-3, 'Y':-1, 'V': 4, 'B':-3, 'Z':-2, 'X':-1, '*':-4},
        'B': {'A':-2, 'R':-1, 'N': 3, 'D': 4, 'C':-3, 'Q': 0, 'E': 1, 'G':-1, 'H': 0, 'I':-3, 'L':-4, 'K': 0,
              'M':-3, 'F':-3, 'P':-2, 'S': 0, 'T':-1, 'W':-4, 'Y':-3, 'V':-3, 'B': 4, 'Z': 1, 'X':-1, '*':-4},
        'Z': {'A':-1, 'R': 0, 'N': 0, 'D': 1, 'C':-3, 'Q': 3, 'E': 4, 'G':-2, 'H': 0, 'I':-3, 'L':-3, 'K': 1,
              'M':-1, 'F':-3, 'P':-1, 'S': 0, 'T':-1, 'W':-3, 'Y':-2, 'V':-2, 'B': 1, 'Z': 4, 'X':-1, '*':-4},
        'X': {'A': 0, 'R':-1, 'N':-1, 'D':-1, 'C':-2, 'Q':-1, 'E':-1, 'G':-1, 'H':-1, 'I':-1, 'L':-1, 'K':-1,
              'M':-1, 'F':-1, 'P':-2, 'S': 0, 'T': 0, 'W':-2, 'Y':-1, 'V':-1, 'B':-1, 'Z':-1, 'X':-1, '*':-4},
        '*': {'A':-4, 'R':-4, 'N':-4, 'D':-4, 'C':-4, 'Q':-4, 'E':-4, 'G':-4, 'H':-4, 'I':-4, 'L':-4, 'K':-4,
              'M':-4, 'F':-4, 'P':-4, 'S':-4, 'T':-4, 'W':-4, 'Y':-4, 'V':-4, 'B':-4, 'Z':-4, 'X':-4, '*': 1}}


BASE2INT = {
    "A": 0,
    "C": 1,
    "G": 2,
    "T": 3
}

INT2BASE = ["A", "C", "G", "T"]



#=============================================================================
# Sequence functions
#

class TranslateError (Exception):
    """Exception raised when a codon cannot be translated correctly.

    Attributes:
        aa: The amino-acid sequence string being reverse-translated.
        dna: The original DNA sequence string.
        a: The amino-acid character that triggered the error.
        codon: The DNA codon that did not match.
    """
    def __init__(self, msg, aa, dna, a, codon):
        """Initialize a TranslateError.

        Args:
            msg: Human-readable error message.
            aa: Amino-acid sequence being processed.
            dna: Original DNA sequence.
            a: The amino-acid character at the point of failure.
            codon: The DNA codon at the point of failure.
        """
        Exception.__init__(self, msg)
        self.aa = aa
        self.dna = dna
        self.a = a
        self.codon = codon



def translate(dna, table=CODON_TABLE):
    """Translate a DNA sequence (with gaps) into an amino-acid sequence.

    Codons containing "N" are translated to "X" (unknown amino acid).
    Gap codons "---" are translated to "-".

    Args:
        dna: DNA string whose length must be a multiple of 3.
        table: Codon-to-amino-acid lookup dict (default CODON_TABLE).

    Returns:
        Amino-acid sequence string.

    Raises:
        AssertionError: If len(dna) is not a multiple of 3.
        KeyError: If a codon is not present in the codon table.
    """

    aa = []

    assert len(dna) % 3 == 0, "dna sequence length is not a multiple of 3"

    for i in range(0, len(dna), 3):
        codon = dna[i:i+3].upper()
        if "N" in codon:
            aa.append("X")     # unkown aa
        else:
            aa.append(table[codon])
    return "".join(aa)


def revtranslate(aa, dna, check=False):
    """Reverse-translate an amino-acid sequence (with gaps) back into DNA.

    The original ungapped DNA sequence must be supplied so that the correct
    codons are restored. Gap characters "-" in aa are expanded to "---" in the
    output.

    Args:
        aa: Amino-acid string (may contain "-" gap characters).
        dna: Original ungapped DNA string used to recover codons.
        check: If True, verify that each codon translates back to the
            expected amino acid (default False).

    Returns:
        DNA string with codons matching the amino-acid sequence.

    Raises:
        TranslateError: If check=True and a codon does not translate to the
            expected amino acid.
    """

    seq = []
    i = 0
    for a in aa:
        if a == "-":
            seq.append("---")
        else:
            codon = dna[i:i+3]
            if check and a != CODON_TABLE.get(codon, "X"):
                raise TranslateError("bad translate", aa, dna, a, codon)
            seq.append(codon)
            i += 3
    return "".join(seq)

_comp = {"A":"T", "C":"G", "G":"C", "T":"A", "N":"N",
         "a":"t", "c":"g", "g":"c", "t":"a", "n":"n",
         "R":"Y", "Y":"R", "S":"W", "W":"S", "K":"M", "M":"K",
         "r":"y", "y":"r", "s":"w", "w":"s", "k":"m", "m":"k",
         "B":"V", "V":"B", "D":"H", "H":"D",
         "b":"v", "v":"b", "d":"h", "h":"d"}

def revcomp(seq):
    """Return the reverse complement of a DNA sequence.

    Handles IUPAC ambiguity codes as well as standard A/C/G/T bases (both
    upper and lower case).

    Args:
        seq: DNA sequence string.

    Returns:
        Reverse-complemented DNA sequence string.
    """

    seq2 = []
    for i in range(len(seq)-1, -1, -1):
        seq2.append(_comp[seq[i]])
    return "".join(seq2)


def gcContent(seq):
    """Compute the GC content fraction of a DNA sequence.

    Args:
        seq: DNA sequence string containing A, C, G, and T characters.

    Returns:
        GC fraction as a float in [0.0, 1.0].
    """
    # Inlined util.histDict: build a frequency dict of characters
    hist = {}
    for c in seq:
        hist[c] = hist.get(c, 0) + 1
    total = hist["A"] + hist["C"] + hist["T"] + hist["G"]

    return (hist["C"] + hist["G"]) / float(total)


#=============================================================================
# Kimura sequence mutation model
#

KIMURA_MATRIX = [
    ['r', 's', 'u', 's'],
    ['s', 'r', 's', 'u'],
    ['u', 's', 'r', 's'],
    ['s', 'u', 's', 'r']
]


def evolveKimuraSeq(seq, time, alpha=1, beta=1):
    """Evolve a DNA sequence under the Kimura two-parameter model.

    Each base is independently substituted according to transition (alpha)
    and transversion (beta) rate parameters over the given evolutionary time.

    Args:
        seq: DNA sequence string (uppercase A/C/G/T only).
        time: Evolutionary branch length (substitutions per site).
        alpha: Transition rate parameter (default 1).
        beta: Transversion rate parameter (default 1).

    Returns:
        Evolved DNA sequence string of the same length as seq.

    Raises:
        AssertionError: If substitution probabilities do not sum to one.
    """
    probs = {
        's': .25 * (1 - math.e**(-4 * beta * time)),
        'u': .25 * (1 + math.e**(-4 * beta * time)
                      - 2*math.e**(-2*(alpha+beta)*time))
    }
    probs['r'] =  1 - 2*probs['s'] - probs['u']

    seq2 = []

    for base in seq:
        cdf = 0
        row = KIMURA_MATRIX[BASE2INT[base]]
        pick = random.random()

        for i in range(4):
            cdf += probs[row[i]]
            if cdf >= pick:
                seq2.append(INT2BASE[i])
                break

    assert len(seq2) == len(seq), "probabilities do not add to one"

    return "".join(seq2)


def evolveKimuraBase(base, time, alpha, beta):
    """Evolve a single DNA base under the Kimura two-parameter model.

    Args:
        base: A single DNA base character (A/C/G/T).
        time: Evolutionary branch length.
        alpha: Transition rate parameter.
        beta: Transversion rate parameter.

    Returns:
        The (possibly substituted) DNA base character.

    Raises:
        AssertionError: If substitution probabilities do not sum to one.
    """
    probs = {
        's': .25 * (1 - math.e**(-4 * beta * time)),
        'u': .25 * (1 + math.e**(-4 * beta * time)
                      - 2*math.e**(-2*(alpha+beta)*time))
    }
    probs['r'] =  1 - 2*probs['s'] - probs['u']

    cdf = 0
    row = KIMURA_MATRIX[BASE2INT[base]]
    pick = random.random()

    for i in range(4):
        cdf += probs[row[i]]
        if cdf >= pick:
            return INT2BASE[i]

    assert False, "probabilities do not add to one"
