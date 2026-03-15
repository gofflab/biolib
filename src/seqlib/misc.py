#!/usr/bin/python
"""Miscellaneous utility functions for sequence analysis, data structures, and pretty printing.

Provides tools for nuID encoding/decoding of nucleotide sequences, dictionary sorting,
pretty-printing of nested data structures, ranking/ordering utilities, and basic string
manipulation functions used across the seqlib package.
"""
import sys


#############
#pygr tools
#############
class Annot:
    """Annotation class for pygr data.

    A lightweight container for genomic annotation records used with the pygr
    genome database library.

    Attributes:
        name: Identifier for the annotation (e.g. gene name or transcript ID).
        chr: Chromosome name (e.g. 'chr1').
        strand: Strand orientation ('+' or '-').
        start: 0-based start coordinate of the annotation.
        end: End coordinate of the annotation.
    """
    def __init__(self,name,chr,strand,start,end):
        """Initialises an Annot instance.

        Args:
            name: Identifier for the annotation.
            chr: Chromosome name.
            strand: Strand orientation ('+' or '-').
            start: 0-based start coordinate.
            end: End coordinate.
        """
        self.name=name
        self.chr=chr
        self.strand=strand
        self.start=start
        self.end=end

##################
#nuID implementation for python
###################
def mreplace(s,chararray=['A','C','G','T','U'],newarray=['0','1','2','3','3']):
    """Replaces multiple characters in a string using paired replacement arrays.

    Iterates over corresponding pairs from chararray and newarray, replacing
    each occurrence of chararray[i] with newarray[i] in sequence.  Defaults
    map the nucleotide alphabet (A, C, G, T, U) to single-digit codes used
    by the nuID encoding scheme.

    Args:
        s: Input string to perform replacements on.
        chararray: List of characters (or substrings) to replace.
        newarray: List of replacement characters (or substrings), paired
            positionally with chararray.

    Returns:
        The modified string after all replacements have been applied.
    """
    for a,b in zip(chararray,newarray):
        s=s.replace(a,b)
    return s

def seq2nuID(seq):
    """Converts a DNA or RNA sequence string into its corresponding nuID.

    The nuID (nucleotide identifier) is a compact, base-64-like encoding of a
    nucleotide sequence that encodes both sequence content and a checksum
    character.  This implementation replaces the standard "_" character in the
    code alphabet with "!" to avoid conflicts with SHRiMP alignment output
    parsing.

    Args:
        seq: A DNA or RNA sequence string (case-insensitive; 'U' is treated
            identically to 'T').

    Returns:
        A nuID string whose first character encodes checksum and padding
        information and whose remaining characters encode successive triplets
        of nucleotides in base-64 space.
    """

    """
        Default code includes "_" as char.  This conflicts with parsing for shrimp.  So for my specific instance,
        "_" has been replaced with "!"
    """
    code = [chr(x) for x in range(65,91)]+[chr(x) for x in range(97,123)]+[str(x) for x in range(0,10)]+[str(x) for x in ("!",".")]
    seq=seq.upper()
    num=mreplace(seq)
    if len(num)%3!=0:
        appLen = 3-len(num)%3
        num = num+(appLen*str(0))
    else:
        appLen=0
    numArray=[]
    checkSumArray=[]
    charArray=[]
    for i in range(0,len(num),3):
        subnum=num[i:i+3]
        numArray.append(subnum)
        code64=int(subnum[0])*4**2+int(subnum[1])*4+int(subnum[2])
        checkSumArray.append(code64)
        charArray.append(code[code64])
    checkSum=sum(checkSumArray)
    res=checkSum%21
    checkCode=code[res*3+appLen]
    #print numarray
    #print checkSumArray
    id = str(checkCode)+"".join(charArray)
    return id

def nuID2seq(nuID):
    """Decodes a nuID string back into the original nucleotide sequence.

    Reverses the nuID encoding produced by seq2nuID.  The first character of
    the nuID encodes checksum and padding length; the remaining characters are
    decoded from base-64 triplets back to the ACGT alphabet.  This
    implementation uses "!" instead of "_" in the code alphabet (matching
    seq2nuID) to avoid conflicts with SHRiMP output parsing.

    Args:
        nuID: A nuID string as produced by seq2nuID.

    Returns:
        The original DNA sequence string (uppercase ACGT).

    Raises:
        AssertionError: If the nuID contains the '.' character as a check code
            (which would indicate a coding error or invalid nuID), or if the
            checksum validation fails.
    """
    """
        Default code includes "_" as char.  This conflicts with parsing for shrimp.  So for my specific instance,
        "_" has been replaced with "!"
    """
    import math
    code = [chr(x) for x in range(65,91)]+[chr(x) for x in range(97,123)]+[str(x) for x in range(0,10)]+[str(x) for x in ("!",".")]
    ind=range(1,len(code)+1)
    names=dict(zip(code,ind))
    numArray=[]
    for l in nuID:
        numArray.append(names[l]-1) #Possibly need to add -1?
    checkCode=int(numArray.pop(0))
    if checkCode==63:
        assert "Coding error or not a nuID!\nCheck code should not include '.'!"
    cutlen=checkCode%3
    res=int(math.floor(checkCode/3))
    num=''
    for i in numArray:
        subDecode = [int(math.floor(i/4**2)),int(math.floor((i%4**2)/4)),int(i%4)]
        newsub = "".join(str(j) for j in subDecode)
        num=num+newsub
    checkSum=sum(numArray)
    if res != checkSum%21:
        assert "Coding Error or not a nuID"
    nucleotide=["A","C","G","T"]
    seq=mreplace(num,['0','1','2','3'],nucleotide)
    seq=seq[:-1]
    return seq

######
#
#Dictionary tools
#
#######

def sort_by_value(d):
    """ Returns the keys of dictionary d sorted by their values """
    items=d.items()
    backitems=[ [v[1],v[0]] for v in items]
    backitems.sort(reverse=True)
    return [ backitems[i][1] for i in range(0,len(backitems))]

def sbv2(d,reverse=False):
    """Returns dictionary items sorted by value, using itemgetter (PEP 265 approach).

    Args:
        d: A dictionary to sort.
        reverse: Not currently used; items are always sorted in descending
            order by value regardless of this parameter.

    Returns:
        A list of (key, value) tuples sorted by value in descending order.
    """
    from operator import itemgetter
    return sorted(d.items(), key=itemgetter(1), reverse=True)

def sortListofDicts(fieldname):
    """useful for sorting a list of dictionaries by a given key (fieldname)
    usage:
    mylist.sort(key=sortListofDicts('start'))  #will sort a list of intervals by i['start']
    """
    return lambda x: x[fieldname]

def sort_dict(d,reverse=True):
    """Returns dictionary items sorted first by value then by key.

    Args:
        d: A dictionary to sort.
        reverse: If True (default), sort in descending order; if False,
            sort in ascending order.

    Returns:
        A list of (key, value) tuples sorted by (value, key) using the
        specified direction.
    """
    return sorted(d.items(), key=lambda item: (item[1], item[0]), reverse=reverse)

########
#
#Pretty Printing
#
########
def pretty_print(f, d, level=-1, maxw=0, maxh=0, gap="", first_gap='', last_gap=''):
    """Recursively pretty-prints a nested Python data structure to a file stream.

    Handles lists, tuples, dicts, class instances, and scalar values, printing
    each with indentation that reflects the nesting depth.  Optionally limits
    the depth of recursion, the width of each printed line, and the number of
    elements printed per container.

    Args:
        f: Output file stream (e.g. sys.stdout or an open file handle).
        d: The data structure to print.
        level: Maximum recursion depth.  -1 (default) means unlimited depth.
            0 means stop recursing and print a repr of the current element.
        maxw: Maximum character width for a single printed line.  0 (default)
            means no width limit.
        maxh: Maximum number of elements to print from any list, tuple, or
            dict at any recursion level.  0 (default) means no limit.
        gap: Indentation prefix inserted before each element inside a
            container.
        first_gap: Prefix printed before the opening bracket/brace/paren of
            a container, or before a scalar value.
        last_gap: Prefix printed before the closing bracket/brace/paren of
            a container.
    """
    # depending on the type of expression, it recurses through its elements
    # and prints with appropriate indentation

    # f   is the output file stream
    # d   is the data structure
    #
    # level is the number of allowed recursive calls, the depth at which
    #       the data structure is explored
    #       default: -1 means never stop recursing early
    # maxw  is the maximum width that will be printed from the last element
    #       of the recursion (when no further recursion is possible, or
    #       the maximal depth has been reached)
    #       default: 0 means every line will be printed in its entirety, regardless
    #                of how long it may be
    # maxh  (max height) is the maximum number of elements that will be
    #       printed from a list or a dictionary, at any level or recursion
    #       default: 0 means every list or dictionary will have all its elements
    #                printed, even if it contains thousands of elements
    #
    # gap is the gap to include before every element of a list/dic/tuple
    # first_gap is the opening gap before the opening bracket, parens or curly braces
    # first_gap is the closing gap before the closing bracket, parens or curly braces

    if level == 0:
        if not isinstance(d, str): d = repr(d)

        if maxw and len(d) > maxw:
            final = ifab(maxw > 20, 10, maxw/2)
            f.write(first_gap+d[:maxw-final]+'...'+d[-final:]+' (%s chars)\n' % len(d))
        else: f.write(first_gap+d+'\n')
    elif isinstance(d, list):
        if not d:
            f.write(first_gap+"[]\n")
            return
        # recurse on lists
        f.write(first_gap+"[\n")
        h = 0
        for el in d:
            pretty_print(f, el, level-1, maxw, maxh, gap+'   ', gap+' ->', gap+'   ')
            if maxh:
                h = h+1
                if h >= maxh and maxh<len(d):
                    f.write(gap+' -> ... (%s in list)\n'%len(d))
                    break
        f.write(last_gap+"]\n")
    elif isinstance(d, tuple):
        if not d:
            f.write(first_gap+"()\n")
            return
        # recurse on tuples
        f.write(first_gap+"(\n")
        h = 0
        for el in d:
            pretty_print(f, el,
                         level     = level-1,
                         maxw      = maxw,
                         maxh      = maxh,
                         gap       = gap+'   ',
                         first_gap = gap+' =>',
                         last_gap  = gap+'   ')
            if maxh:
                h = h+1
                if h >= maxh and maxh<len(d):
                    f.write(gap+' => ... (%s in tuple)\n'%len(d))
                    break
        f.write(last_gap+")\n")
    elif isinstance(d, dict):
        if not d:
            f.write(first_gap+"{}\n")
            return
        # recurse on dictionaries
        f.write(first_gap+"{\n")
        keys = sorted(d.keys())
        key_strings = [ifab(isinstance(k, str), k, repr(k)) for k in keys]
        maxlen = max(map(len, key_strings))
        h = 0
        for k,key_string in zip(keys, key_strings):
            key_string = sfill(key_string,maxlen,'.')
            blank_string = ' '*len(key_string)
            pretty_print(f, d[k],
                         level     = level-1,
                         maxw      = maxw,
                         maxh      = maxh,
                         gap       = gap+'    %s'%blank_string,
                         first_gap = gap+'  %s: '%key_string,
                         last_gap  = gap+'    %s'%blank_string)
            if maxh:
                h = h+1
                if h >= maxh and maxh<len(keys):
                    remaining_keys = []
                    for k in keys[h:]:
                        if isinstance(k, tuple):
                            remaining_keys.append(repr(k))
                        else:
                            remaining_keys.append('%s'%k)
                    remaining_keys = ','.join(remaining_keys)
                    #f.write(gap+'  %s (%s keys)\n'%(remaining_keys, len(keys)))
                    pretty_print(f, '  %s (%s keys)'%(remaining_keys, len(keys)),0,maxw,0,
                                 gap,gap,'')
                    break

            #gap+' '*(len(key_string)+3), '', gap+' '*(len(key_string)+5))
        f.write(last_gap+"}\n")
    elif hasattr(d, '__dict__') and not isinstance(d, (list, tuple, dict, str, int, float, bool)):
        fields = dir(d)

        if not fields:
            f.write(first_gap+"*EmptyClass*\n")
            return
        # recurse on classes
        f.write(first_gap+"*ClassInstance %s\n"%d)
        fields.sort()
        key_strings = [ifab(isinstance(k, str), k, repr(k)) for k in fields]
        maxlen = max(map(len, key_strings))
        h = 0
        for k,key_string in zip(fields, key_strings):
            key_string = sfill(key_string,maxlen,'.')
            blank_string = ' '*len(key_string)
            pretty_print(f, eval('d.'+k),
                         level     = level-1,
                         maxw      = maxw,
                         maxh      = maxh,
                         gap       = gap+'    %s'%blank_string,
                         first_gap = gap+'  %s: '%key_string,
                         last_gap  = gap+'    %s'%blank_string)
            if maxh:
                h = h+1
                if h >= maxh and maxh<len(keys):
                    remaining_keys = []
                    for k in keys[h:]:
                        if isinstance(k, tuple):
                            remaining_keys.append(repr(k))
                        else:
                            remaining_keys.append('%s'%k)
                    remaining_keys = ','.join(remaining_keys)
                    #f.write(gap+'  %s (%s keys)\n'%(remaining_keys, len(keys)))
                    pretty_print(f,
                                 '  %s (%s keys)'%(remaining_keys, len(keys)),
                                 0,
                                 maxw,
                                 0,
                                 gap,
                                 gap,
                                 '')
                    break

            #gap+' '*(len(key_string)+3), '', gap+' '*(len(key_string)+5))
        f.write(last_gap+"*\n")
    elif type(d) == type(""):
        # simply print strings (no quotes)
        if maxw and len(d)>maxw:
            final = ifab(maxw > 20, 10, maxw/2)
            f.write(first_gap+d[:maxw-final]+'..'+d[-final:]+' (%s)\n' % len(d))
        else:
            f.write(first_gap+d+'\n')
    else:
        # string conversion of all other types
        if maxw and len(repr(d))>maxw:
            final = ifab(maxw > 20, 10, maxw/2)
            f.write(first_gap+repr(d)[:maxw-final]+'..'+repr(d)[-final:]+' (%s)\n' % len(repr(d)))
        else:
            f.write(first_gap+repr(d)+'\n')

def pp(d,level=-1,maxw=0,maxh=0,parsable=0):
    """Pretty-prints a data structure to stdout.

    Wrapper around pretty_print that writes to sys.stdout.  When parsable is
    set to a truthy value the standard library pprint module is used instead,
    which produces output that can be eval'd back to the original structure.

    Args:
        d: The data structure to print.
        level: Maximum recursion depth passed to pretty_print.  -1 means
            unlimited.
        maxw: Maximum line width passed to pretty_print (or pprint width when
            parsable is set).  0 means no limit.
        maxh: Maximum container height passed to pretty_print.  0 means no
            limit.
        parsable: If 0 (default), use pretty_print for human-readable output.
            If non-zero, use the standard library pprint module.
    """
    if not parsable:
        pretty_print(sys.stdout, d, level, maxw, maxh, '', '', '')
    else:
        import pprint
        if maxw: pp2 = pprint.PrettyPrinter(width=maxw, indent=1)#, depth=level
        else: pp2 = pprint.PrettyPrinter(indent=1)#, depth=level
        pp2.pprint(d)

def test_pp():
    """Runs a self-contained smoke test of the pp / pretty_print functions.

    Calls pp with a heterogeneous nested data structure containing dicts,
    lists, tuples, integers, strings, and a lambda.  Output is written to
    stdout.  No return value.
    """
    pp({'one': ('two',3,[4,5,6]),
        7: (lambda x: 8*9),
        'ten': ['ele', {'ven': 12,
                        (13,14): '15'}]})

###################################
#
#Boolean Functions
#
####################################
def ifab(test, a, b):
    """x = ifab(test, a, b)
       WARNING:  Both 'a' and 'b' are evaluated
       C equivalent: x = test?a:b;
       Scheme equiv: (set x (if test a b))
       Python equiv: test and a or b
       None of the equivalents evaluates both arguments
    """
    if test: return a
    else: return b


###################################
#
#String Functions
#
####################################
def sfill(s, length, fill_char = '.'):
    """Pads a string on the right with a fill character until it reaches the target length.

    Example::

        sfill('hello', 18, '.') -> 'hello.............'
        #                           <---  18 chars  --->

    Useful for aligning dictionary keys when pretty-printing:
    ``one......: 1``, ``five.....: 5``, ``seventeen: 17``.

    Args:
        s: The input string to pad.
        length: The desired total length of the returned string.
        fill_char: The character used for padding (default: '.').

    Returns:
        The input string right-padded with fill_char to the specified length.
        If the input string is already at least as long as length, it is
        returned unchanged.
    """
    #  Appends fill_char to the string s until it reaches length length
    #  ex:  sfill('hello',18,'.') -> hello...............
    #                                <---  18 chars  --->
    # useful for printing dictionaries in a cute way
    #    one......: 1
    #    five.....: 5
    #    seventeen: 17


    #list = map(None, s)
    #list.extend(map(None, fill_char*(length - len(list))))
    #return string.join(list, '')

    return s + fill_char*(length-len(s))

def rstrips(s, suffix):
    """Strips a specific suffix from the right end of a string.

    Unlike str.rstrip, this function removes the exact suffix string rather
    than a set of characters.

    Args:
        s: The input string.
        suffix: The exact suffix to remove.  If empty or not present at the
            end of s, the string is returned unchanged.

    Returns:
        The input string with the suffix removed from the right end, or the
        original string if the suffix was not found.
    """
    if suffix and s.endswith(suffix):
        s = s[:-len(suffix)]
    return s

def hamming_distance(s1, s2):
    """Returns the hamming (or edit) distance between two strings or list of iterable elements"""
    assert len(s1) == len(s2)
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

######################################
#
#Ranking and Ordering
#
######################################
from random import sample  # noqa: E402


def order(x, NoneIsLast = True, decreasing = False):
    """
    Returns the ordering of the elements of x. The list
    [ x[j] for j in order(x) ] is a sorted version of x.

    Missing values in x are indicated by None. If NoneIsLast is true,
    then missing values are ordered to be at the end.
    Otherwise, they are ordered at the beginning.
    """
    omitNone = False
    if NoneIsLast == None:
        NoneIsLast = True
        omitNone = True

    n  = len(x)
    ix = range(n)
    if None not in x:
        ix.sort(reverse = decreasing, key = lambda j : x[j])
    else:
        # Handle None values properly.
        def key(i, x = x):
            elem = x[i]
            # Valid values are True or False only.
            if decreasing == NoneIsLast:
                return elem is not None, elem
            else:
                return elem is None, elem
        ix = range(n)
        ix.sort(key=key, reverse=decreasing)

    if omitNone:
        n = len(x)
        for i in range(n-1, -1, -1):
            if x[ix[i]] == None:
                n -= 1
        return ix[:n]
    return ix


def rank(x, NoneIsLast=True, decreasing = False, ties = "first"):
    """
    Returns the ranking of the elements of x. The position of the first
    element in the original vector is rank[0] in the sorted vector.

    Missing values are indicated by None.  Calls the order() function.
    Ties are NOT averaged by default. Choices are:
         "first" "average" "min" "max" "random" "average"
    """
    omitNone = False
    if NoneIsLast == None:
        NoneIsLast = True
        omitNone = True
    O = order(x, NoneIsLast = NoneIsLast, decreasing = decreasing)
    R = O[:]
    n = len(O)
    for i in range(n):
        R[O[i]] = i
    if ties == "first" or ties not in ["first", "average", "min", "max", "random"]:
        return R

    blocks     = []
    isnewblock = True
    newblock   = []
    for i in range(1,n) :
        if x[O[i]] == x[O[i-1]]:
            if i-1 not in newblock:
                newblock.append(i-1)
            newblock.append(i)
        else:
            if len(newblock) > 0:
                blocks.append(newblock)
                newblock = []
    if len(newblock) > 0:
        blocks.append(newblock)

    for i, block  in enumerate(blocks):
        # Don't process blocks of None values.
        if x[O[block[0]]] == None:
            continue
        if ties == "average":
            s = 0.0
            for j in block:
                s += j
            s /= float(len(block))
            for j in block:
                R[O[j]] = s
        elif ties == "min":
            s = min(block)
            for j in block:
                R[O[j]] = s
        elif ties == "max":
            s =max(block)
            for j in block:
                R[O[j]] = s
        elif ties == "random":
            s = sample([O[i] for i in block], len(block))
            for i,j in enumerate(block):
                R[O[j]] = s[i]
        else:
            for i,j in enumerate(block):
                R[O[j]] = j
    if omitNone:
        R = [ R[j] for j in range(n) if x[j] != None]
    return R

def uniqify(seq):
    """Returns the unique elements of an iterable as a list.

    Not order-preserving: the returned list may appear in arbitrary order
    because uniqueness is tracked via a dictionary.

    Args:
        seq: An iterable of hashable elements.

    Returns:
        A list containing each unique element from seq exactly once.
    """
    # Not order preserving
    keys = {}
    for e in seq:
        keys[e] = 1
    return list(keys.keys())
