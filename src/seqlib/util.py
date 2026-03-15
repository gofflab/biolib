"""

    Common Utilities

    file: util.py 
    authors: Matt Rasmussen
    date: 11/30/05
 
    Provides basic functional programming functions for manipulating lists and 
    dicts.  Also provides common utilities (timers, plotting, histograms)
    
"""



# python libs
import copy
import math
import os
import re
import sys
from functools import cmp_to_key

#
# see bottom of file for other imports
#


# Note: I had trouble using 1e1000 directly, because bytecode had trouble
# representing infinity (possibly)
INF = float("1e1000")


# Python 3 compatibility: cmp() was removed
def cmp(a, b):
    """Three-way comparison function for Python 3 compatibility.

    Args:
        a: First value.
        b: Second value.

    Returns:
        1 if a > b, -1 if a < b, 0 if equal.
    """
    return (a > b) - (a < b)




class Bundle (dict):
    """A small class for creating a closure of variables.

    Handy for nested functions that need to assign to variables in an outer
    scope. Attributes and dictionary keys are kept in sync.

    Example::

        def func1():
            this = Bundle(var1=0, var2="hello")
            def func2():
                this.var1 += 1
            func2()
            print(this.var1)
        func1()
        # prints: 1
    """

    def __init__(self, **variables):
        """Initialize a Bundle with keyword arguments as attributes.

        Args:
            **variables: Arbitrary keyword arguments that become both
                attributes (self.key) and dictionary entries.
        """
        for key, val in variables.items():
            setattr(self, key, val)
            dict.__setitem__(self, key, val)

    def __setitem__(self, key, val):
        """Set a key both as an attribute and as a dict entry.

        Args:
            key: Attribute/key name.
            val: Value to assign.
        """
        setattr(self, key, val)
        dict.__setitem__(self, key, val)



class Dict (dict):
    """A nested dictionary with configurable dimensionality and default values.

    Accessing a missing key returns (and optionally inserts) a default value
    or a nested Dict of one lower dimension, enabling multi-dimensional sparse
    containers without explicit initialisation.
    """


    def __init__(self, items=None, dim=1, default=None, insert=True):
        """Initialize a Dict.

        Args:
            items: Initial items to populate the dict (dict, list of pairs,
                or other iterable). If an int is passed, it is treated as
                the old-style positional dim argument for backwards
                compatibility.
            dim: Number of nesting dimensions (default 1).
            default: Default value returned for missing leaf-level keys
                (default None).
            insert: If True, accessing a missing key inserts the default
                value automatically (default True).
        """

        if isinstance(items, int):
            # backwards compatiability
            default = dim
            dim = items
        elif items is not None:
            dict.__init__(self, items)

        self._dim = dim
        self._null = default
        self._insert = insert

        # backwards compatiability
        self.data = self


    def __getitem__(self, i):
        """Return the value for key i, inserting a default if missing.

        Args:
            i: The key to look up.

        Returns:
            The stored value, or a default Dict/copy of null if the key was
            absent.
        """
        if i not in self:
            if self._dim > 1:
                ret = Dict(self._dim - 1, self._null)
            else:
                ret = copy.copy(self._null)
            if self._insert:
                self[i] = ret
            return ret
        return dict.__getitem__(self, i)


    def has_keys(self, *keys):
        """Check whether a sequence of nested keys all exist.

        Args:
            *keys: Keys to check at successive nesting levels.

        Returns:
            True if all keys are present at the corresponding nesting levels.
        """
        if len(keys) == 0:
            return True
        elif len(keys) == 1:
            return keys[0] in self
        else:
            return keys[0] in self and \
                   self[keys[0]].has_keys(*keys[1:])

    def write(self, out = sys.stdout):
        """Write a human-readable representation of the dict to a stream.

        Args:
            out: Output stream to write to (default sys.stdout).
        """
        def walk(node, path):
            if node.dim == 1:
                for i in node:
                    out.write("  ")
                    for j in path:
                        out.write(str(j) + ", ")
                    print(i, ":", node[i], file=out)
            else:
                for i in node:
                    walk(node[i], path + [i])

        print("< DictMatrix", file=out)
        walk(self, [])
        print(">", file=out)




class Percent (float):
    """A float subclass that formats itself as a percentage string.

    Attributes:
        digits: Number of decimal places used when formatting (default 1).
    """
    digits = 1

    def __str__(self):
        """Return the value formatted as a percentage with self.digits decimals.

        Returns:
            String such as "42.0" representing 42.0% (i.e. float value 0.42).
        """
        return (("%%.%df" % self.digits) % (float(self) * 100))

    def __repr__(self):
        """Return the same string as __str__."""
        return str(self)


class PushIter (object):
    """An iterator wrapper that allows pushing items back to the front of the stream.

    Wraps any iterable and provides a push() method to prepend items.
    """

    def __init__(self, it):
        """Initialize a PushIter from any iterable.

        Args:
            it: Any iterable to wrap.
        """
        self._it = iter(it)
        self._queue = []

    def __iter__(self):
        """Return self as the iterator."""
        return self

    def __next__(self):
        """Return the next item, preferring items from the push queue.

        Returns:
            The next item from the queue if non-empty, otherwise from the
            underlying iterator.
        """
        if len(self._queue) > 0:
            return self._queue.pop()
        else:
            return self.next(_it)

    def push(self, item):
        """Push a new item onto the front of the iteration stream.

        Args:
            item: Item to prepend to the iteration.
        """
        self._queue.append(item)


def exceptDefault(func, val, exc=Exception):
    """Call func() and return val if the specified exception is raised.

    Args:
        func: A zero-argument callable to invoke.
        val: Default value to return on exception.
        exc: Exception type (or tuple of types) to catch (default Exception).

    Returns:
        The return value of func(), or val if exc was raised.
    """
    try:
        return func()
    except exc:
        return val


#=============================================================================
# list and dict functions for functional programming

def equal(* vals):
    """Returns True if all arguments are equal"""
    if len(vals) < 2:
        return True
    a = vals[0]
    for b in vals[1:]:
        if a != b:
            return False
    return True


def remove(lst, *vals):
    """Returns a copy of list 'lst' with values 'vals' removed
    """
    lst2 = []
    delset = set(vals)
    for i in lst:
        if i not in delset:
            lst2.append(i)
    return lst2


def sort(lst, compare=None, key=None, reverse=False):
    """Returns a sorted copy of a list
       
       lst     -- a list to sort
       compare -- a comparison function (deprecated in Python 3, use key=)
       key     -- function of one arg to map items
       reverse -- when True reverse sorting
    """
    lst2 = list(lst)
    if compare is not None and compare is not cmp:
        lst2.sort(key=cmp_to_key(compare), reverse=reverse)
    else:
        lst2.sort(key=key, reverse=reverse)
    return lst2


def reverse(lst):
    """Returns a reversed copy of a list
    """
    lst2 = list(lst)
    lst2.reverse()
    return lst2


def cget(mat, *i):
    """Returns the column(s) '*i' of a 2D list 'mat'
        
       mat -- matrix or 2D list 
       *i  -- columns to extract from matrix
       
       notes:
       If one column is given, the column is returned as a list.
       If multiple columns are given, a list of columns (also lists) is returned
    """

    if len(i) == 1:
        return [row[i[0]] for row in mat]
    else:
        return [[row[index] for row in mat]
                for index in i]


def mget(lst, ind):
    """Returns a list 'lst2' such that lst2[i] = lst[ind[i]]
       
       Or in otherwords, get the subsequence 'lst'       
    """
    return [lst[i] for i in ind]



def concat(* lists):
    """Concatenates several lists into one
    """

    lst = []
    for l in lists:
        lst.extend(l)
    return lst



def subdict(dic, keys):
    """
    Returns a new dictionary dic2 such that
    dic2[i] = dic[i] for all i in keys

    dic  -- a dictionary
    keys -- a list of keys
    """
    dic2 = {}
    for key in keys:
        if key in dic:
            dic2[key] = dic[key]
    return dic2


def revdict(dic, allowdups=False):
    """
    Reverses a dict 'dic' such that the keys become values and the 
    values become keys.
    
    allowdups -- if True, one of several key-value pairs with the same value 
                 will be arbitrarily choosen.  Otherwise an expection is raised
    """

    dic2 = {}
    if allowdups:
        for key, val in dic.items():
            dic2[val] = key
    else:
        for key, val in dic.items():
            assert key not in dic2, "duplicate value '%s' in dict" % val
            dic2[val] = key

    return dic2


def list2lookup(lst):
    """
    Creates a dict where each key is lst[i] and value is i
    """

    lookup = {}
    for i in range(len(lst)):
        lookup[lst[i]] = i
    return lookup


def mapdict(dic, key=lambda x: x, val=lambda x: x,
            keyfunc=None, valfunc=None):
    """
    Creates a new dict where keys and values are mapped
    
    keyfunc and valfunc are DEPRECATED
    
    """

    if keyfunc is not None:
        key = keyfunc
    if valfunc is not None:
        val = valfunc

    dic2 = {}
    for k, v in dic.items():
        dic2[key(k)] = val(v)

    return dic2


def mapwindow(func, size, lst):
    """Apply a function 'func' to a sliding window of size 'size' within
       a list 'lst'"""
    lst2 = []
    lstlen = len(lst)
    radius = int(size // 2)

    for i in range(lstlen):
        radius2 = min(i, lstlen - i - 1, radius)
        lst2.append(func(lst[i-radius2:i+radius2+1]))

    return lst2


def groupby(func, lst, multi=False):
    """Places i and j of 'lst' into the same group if func(i) == func(j).
       
       func -- is a function of one argument that maps items to group objects
       lst  -- is a list of items
       multi -- if True, func must return a list of keys (key1, ..., keyn) for
                item a.  groupby will retrun a nested dict 'dct' such that
                dct[key1]...[keyn] == a
       
       returns:
       a dictionary such that the keys are groups and values are items found in
       that group
    """

    if not multi:
        dct = {}
        for i in lst:
            dct.setdefault(func(i), []).append(i)
    else:
        dct = {}
        for i in lst:
            keys = func(i)
            d = dct
            for key in keys[:-1]:
                d = d.setdefault(key, {})
            d.setdefault(keys[-1], []).append(i)

    return dct


def unique(lst):
    """
    Returns a copy of 'lst' with only unique entries.
    The list is stable (the first occurance is kept).
    """

    found = set()

    lst2 = []
    for i in lst:
        if i not in found:
            lst2.append(i)
            found.add(i)

    return lst2


def flatten(lst, depth=INF):
    """
    Flattens nested lists/tuples into one list
    
    depth -- specifies how deep flattening should occur
    """

    flat = []

    for elm in lst:
        if hasattr(elm, "__iter__") and depth > 0:
            flat.extend(flatten(elm, depth-1))
        else:
            flat.append(elm)

    return flat


def mapapply(funcs, lst):
    """
    apply each function in 'funcs' to one element in 'lst'
    """

    lst2 = []
    for func, item in zip(funcs, lst):
        lst2.append(func(item))
    return lst2


def cumsum(vals):
    """Returns a cumalative sum of vals (as a list)"""

    lst = []
    tot = 0
    for v in vals:
        tot += v
        lst.append(tot)
    return lst

def icumsum(vals):
    """Returns a cumalative sum of vals (as an iterator)"""

    tot = 0
    for v in vals:
        tot += v
        yield tot


def frange(start, end, step):
    """Generates a range of floats 
    
       start -- begining of range
       end   -- end of range
       step  -- step size
    """

    i = 0
    val = start
    while val < end:
        yield val
        i += 1
        val = start + i * step





#=============================================================================
# simple matrix functions

def make_matrix(nrows, ncols, val = 0):
    """Create a 2D list (matrix) with given dimensions and a fill value.

    Args:
        nrows: Number of rows.
        ncols: Number of columns.
        val: Fill value for each cell (default 0); each cell gets a copy.

    Returns:
        A list of lists of shape (nrows, ncols) filled with copies of val.
    """
    mat = []
    for i in range(nrows):
        row = []
        mat.append(row)
        for j in range(ncols):
            row.append(copy.copy(val))
    return mat
makeMatrix = make_matrix


def transpose(mat):
    """
    Transpose a matrix
    
    Works better than zip() in that rows are lists not tuples
    """

    assert equal(* map(len, mat)), "rows are not equal length"

    mat2 = []

    for j in range(len(mat[0])):
        row2 = []
        mat2.append(row2)
        for row in mat:
            row2.append(row[j])

    return mat2


def submatrix(mat, rows=None, cols=None):
    """
    Returns a submatrix of 'mat' with only the rows and columns specified

    Rows and columns will appear in the order as indicated in 'rows' and 'cols'
    """

    if rows == None:
        rows = range(len(mat))
    if cols == None:
        cols = range(len(mat[0]))

    mat2 = []

    for i in rows:
        newrow = []
        mat2.append(newrow)
        for j in cols:
            newrow.append(mat[i][j])

    return mat2


def map2(func, *matrix):
    """
    Maps a function onto the elements of a matrix
    
    Also accepts multiple matrices.  Thus matrix addition is
    
    map2(add, matrix1, matrix2)
    
    """

    matrix2 = []

    for i in range(len(matrix[0])):
        row2 = []
        matrix2.append(row2)

        for j in range(len(matrix[0][i])):
            args = [x[i][j] for x in matrix]
            row2.append(func(* args))

    return matrix2


def min2(matrix):
    """Finds the minimum of a 2D list or matrix
    """
    return min(map(min, matrix))


def max2(matrix):
    """Finds the maximum of a 2D list or matrix
    """
    return max(map(max, matrix))


def range2(width, height):
    """Iterates over the indices of a matrix
    
       Thus list(range2(3, 2)) returns
        [(0, 0), (0, 1), (1, 0), (1, 1), (2, 0), (2, 1)]
    """

    for i in range(width):
        for j in range(height):
            yield i, j


#=============================================================================
# list counting and finding functions


def count(func, lst):
    """
    Counts the number of times func(x) is True for x in list 'lst'
       
    See also:
        counteq(a, lst)   count items equal to a
        countneq(a, lst)  count items not equal to a
        countle(a, lst)   count items less than or equal to a
        countlt(a, lst)   count items less than a
        countge(a, lst)   count items greater than or equal to a
        countgt(a, lst)   count items greater than a
    """
    n = 0
    for i in lst:
        if func(i):
            n += 1
    return n

def counteq(a, lst):
    """Count items in lst equal to a."""
    return count(eqfunc(a), lst)

def countneq(a, lst):
    """Count items in lst not equal to a."""
    return count(neqfunc(a), lst)

def countle(a, lst):
    """Count items in lst less than or equal to a."""
    return count(lefunc(a), lst)

def countlt(a, lst):
    """Count items in lst strictly less than a."""
    return count(ltfunc(a), lst)

def countge(a, lst):
    """Count items in lst greater than or equal to a."""
    return count(gefunc(a), lst)

def countgt(a, lst):
    """Count items in lst strictly greater than a."""
    return count(gtfunc(a), lst)


def find(func, *lsts):
    """
    Returns the indices 'i' of 'lst' where func(lst[i]) == True
    
    if N lists are passed, N arguments are passed to 'func' at a time.
    Thus, find(func, list1, list2) returns the list of indices 'i' where 
    func(list1[i], list2[i]) == True
    
    See also:
        findeq(a, lst)   find items equal to a
        findneq(a, lst)  find items not equal to a
        findle(a, lst)   find items less than or equal to a
        findlt(a, lst)   find items less than a
        findge(a, lst)   find items greater than or equal to a
        findgt(a, lst)   find items greater than a
    """

    pos = []

    if len(lsts) == 1:
        # simple case, one list
        lst = lsts[0]
        for i in range(len(lst)):
            if func(lst[i]):
                pos.append(i)
    else:
        # multiple lists given
        assert equal(* map(len, lsts)), "lists are not same length"

        #nvars = len(lsts)
        for i in range(len(lsts[0])):
            if func(* [x[i] for x in lsts]):
                pos.append(i)

    return pos

def findeq(a, lst):
    """Return indices of items in lst equal to a."""
    return find(eqfunc(a), lst)

def findneq(a, lst):
    """Return indices of items in lst not equal to a."""
    return find(neqfunc(a), lst)

def findle(a, lst):
    """Return indices of items in lst less than or equal to a."""
    return find(lefunc(a), lst)

def findlt(a, lst):
    """Return indices of items in lst strictly less than a."""
    return find(ltfunc(a), lst)

def findge(a, lst):
    """Return indices of items in lst greater than or equal to a."""
    return find(gefunc(a), lst)

def findgt(a, lst):
    """Return indices of items in lst strictly greater than a."""
    return find(gtfunc(a), lst)


def islands(lst):
    """Takes a iterable and returns islands of equal consequtive items 
    
    Return value is a dict with the following format
    
    counts = {elm1: [(start,end), (start,end), ...],
              elm2: [(start,end), (start,end), ...]
              ...}
    
    where for each (start,end) in counts[elm1] we have lst[start:end] only 
    containing elm1
    
    """

    counts = {}
    NULL = Bundle() # unique NULL
    last = NULL
    start = 0

    for i, x in enumerate(lst):
        if x != last and last != NULL:
            counts.setdefault(last, []).append((start, i))
            start = i
        last = x
    if last != NULL:
        counts.setdefault(last, []).append((start, i+1))

    return counts



#=============================================================================
# max and min functions

def argmax(lst, key=lambda x: x):
    """
    Find the index 'i' in 'lst' with maximum lst[i]
    
    lst -- list to search
    key -- function to apply to each lst[i].
           argmax(lst, key=func) --> argmax(map(key, lst))
    """

    assert len(lst) > 0
    top = 0
    topval = key(lst[0])
    for i in range(1, len(lst)):
        val = key(lst[i])
        if val > topval:
            top = i
            topval = val
    return top


def argmin(lst, key=lambda x: x):
    """
    Find the index 'i' in 'lst' with minimum lst[i]
    
    lst -- list to search
    key -- function to apply to each lst[i].
           argmin(lst, key=func) --> argmin(map(key, lst))
    """

    assert len(lst) > 0
    low = 0
    lowval = key(lst[0])
    for i in range(1, len(lst)):
        val = key(lst[i])
        if val < lowval:
            low = i
            lowval = val
    return low

def maxfunc(func, lst):
    """Find the element 'e' in 'lst' with maximum func(e)"""
    top = -INF
    topi = None
    for i in lst:
        val = func(i)
        if val > top:
            top = val
            topi = i
    return topi


def minfunc(func, lst):
    """Find the element 'e' in 'lst' with minimum func(e)"""
    low = INF
    lowi = None
    for i in lst:
        val = func(i)
        if val < low:
            low = val
            lowi = i
    return lowi



#=============================================================================
# math functions

#
# comparison function factories
#
# These functions will return convenient comparison functions.
#
# example:
#   filter(ltfunc(4), lst) ==> returns all values in lst less than 4
#   count(ltfunc(4), lst)  ==> returns the number of values in lst < 4
#

def eqfunc(a):
    """Return a function that tests equality with a."""
    return lambda x: x == a

def neqfunc(a):
    """Return a function that tests inequality with a."""
    return lambda x: x != a

def ltfunc(a):
    """Return a function that tests x < a."""
    return lambda x: x < a

def gtfunc(a):
    """Return a function that tests x > a."""
    return lambda x: x > a

def lefunc(a):
    """Return a function that tests x <= a."""
    return lambda x: x <= a

def gefunc(a):
    """Return a function that tests x >= a."""
    return lambda x: x >= a

def withinfunc(a, b, ainc=True, binc=True):
    """Return a function that tests whether x is within the range [a, b].

    Args:
        a: Lower bound.
        b: Upper bound.
        ainc: If True, the lower bound is inclusive (default True).
        binc: If True, the upper bound is inclusive (default True).

    Returns:
        A one-argument function returning True if x is in the specified range.
    """
    if ainc:
        if binc:
            return lambda x: a <= x <= b
        else:
            return lambda x: a <= x < b
    else:
        if binc:
            return lambda x: a < x <= b
        else:
            return lambda x: a < x < b


def sign(num):
    """Returns the sign of a number"""
    return (num > 0) - (num < 0)

def lg(num):
    """Retruns the log_2 of a number"""
    return math.log(num, 2)

def add(a, b):
    """Return a + b."""
    return a + b

def sub(a, b):
    """Return a - b."""
    return a - b

def mul(a, b):
    """Return a * b."""
    return a * b

def idiv(a, b):
    """Return a / b (true division)."""
    return a / b

def div(a, b):
    """Return a / float(b)."""
    return a / float(b)

def safediv(a, b, default=INF):
    """Divide a by b, returning default on ZeroDivisionError.

    Args:
        a: Numerator.
        b: Denominator.
        default: Value to return when b is zero (default INF).

    Returns:
        a / float(b), or default if b is zero.
    """
    try:
        return a / float(b)
    except ZeroDivisionError:
        return default

def safelog(x, base=math.e, default=-INF):
    """Compute log(x) in the given base, returning default on error.

    Args:
        x: Value to take the logarithm of.
        base: Logarithm base (default math.e for natural log).
        default: Value to return when x <= 0 or overflow occurs (default -INF).

    Returns:
        math.log(x, base), or default on OverflowError or ValueError.
    """
    try:
        return math.log(x, base)
    except (OverflowError, ValueError):
        return default

def invcmp(a, b):
    """Return the reversed comparison of a and b (i.e. cmp(b, a)).

    Args:
        a: First value.
        b: Second value.

    Returns:
        1 if b > a, -1 if b < a, 0 if equal.
    """
    return cmp(b, a)  # cmp is defined locally above

def clamp(x, low, high):
    """Clamps a value 'x' between the values 'low' and 'high'
       If low == None, then there is no lower bound
       If high == None, then there is no upper bound
    """

    if high != None and x > high:
        return high
    elif low != None and x < low:
        return low
    else:
        return x

def clampfunc(low, high):
    """Return a function that clamps its argument between low and high.

    Args:
        low: Lower bound (or None for no lower bound).
        high: Upper bound (or None for no upper bound).

    Returns:
        A one-argument function equivalent to clamp(x, low, high).
    """
    return lambda x: clamp(x, low, high)



def compose2(f, g):
    """
    Compose two functions into one

    compose2(f, g)(x) <==> f(g(x))
    """
    return lambda *args, **kargs: f(g(*args, **kargs))


def compose(*funcs):
    """Composes two or more functions into one function
    
       example:
       compose(f,g,h,i)(x) <==> f(g(h(i(x))))
    """

    funcs = reversed(funcs)
    f = next(funcs)
    for g in funcs:
        f = compose2(g, f)
    return f


def overlap(a, b, x, y, inc=True):
    """
    Returns True if range [a,b] overlaps [x,y]
    
    inc -- if True, treat [a,b] and [x,y] as inclusive
    """
    if inc:
        return (y >= a) and (x <= b)
    else:
        return (y > a) and (x < b)



#=============================================================================
# regex
#

def match(pattern, text):
    """
    A quick way to do pattern matching.
    
    remember: to name tokens use (?P<name>pattern)
    """

    m = re.match(pattern, text)

    if m == None:
        return {}
    else:
        return m.groupdict()


def evalstr(text):
    """Replace expressions in a string (aka string interpolation)

    ex:
    >>> name = 'Matt'
    >>> evalstr("My name is ${name} and my age is ${12+12}")
    'My name is Matt and my age is 24'
    
    "${!expr}" expands to "${expr}"
    
    """

    # get environment of caller
    frame = sys._getframe(1)
    global_dict = frame.f_globals
    local_dict = frame.f_locals

    # find all expression to replace
    m = re.finditer(r"\$\{(?P<expr>[^\}]*)\}", text)

    # build new string
    try:
        strs = []
        last = 0
        for x in m:
            expr = x.groupdict()['expr']

            strs.append(text[last:x.start()])

            if expr.startswith("!"):
                strs.append("${" + expr[1:] + "}")
            else:
                strs.append(str(eval(expr, global_dict, local_dict)))
            last = x.end()
        strs.append(text[last:len(text)])
    except Exception as e:
        raise Exception("evalstr: " + str(e))

    return "".join(strs)


#=============================================================================
# common Input/Output

def read_ints(filename):
    """Read a list of integers from a file (one int per line)
    
       filename may also be a stream
    """

    infile = open_stream(filename)
    vec = []
    for line in infile:
        vec.append(int(line))
    return vec
readInts = read_ints


def read_floats(filename):
    """Read a list of floats from a file (one float per line)
    
       filename may also be a stream
    """
    infile = open_stream(filename)
    vec = []
    for line in infile:
        vec.append(float(line))
    return vec
readFloats = read_floats


def read_strings(filename):
    """Read a list of strings from a file (one string per line)
    
       filename may also be a stream
    """
    infile = open_stream(filename)
    vec = [line.rstrip() for line in infile]
    return vec
readStrings = read_strings

def read_dict(filename, delim="\t", keytype=str, valtype=str):
    """Read a dict from a file
       
       filename may also be a stream
    """

    infile = open_stream(filename)
    dct = {}

    for line in infile:
        tokens = line.rstrip("\n").split(delim)
        assert len(tokens) >= 2, line
        dct[keytype(tokens[0])] = valtype(tokens[1])

    return dct
readDict = read_dict

def write_list(filename, lst):
    """Write a list of anything (ints, floats, strings, etc) to a file.
    
       filename may also be a stream
    """
    out = open_stream(filename, "w")
    for i in lst:
        print(i, file=out)
writeList = write_list
writeVector = write_list


def write_dict(filename, dct, delim="\t"):
    """Write a dictionary to a file"""

    out = open_stream(filename, "w")
    for k, v in dct.items():
        out.write("%s%s%s\n" % (str(k), delim, str(v)))
writeDict = write_dict


'''
def makeReopenStream(stream):
    """Object used to wrap a stream that is 'opened' multiple times.
       Will ignore first close"""
    
    closecount = [0]
    
    # save old close
    old_close = stream.close
    
    def new_close():
        closecount[0] += 1
        if closecount[0] > 1:
            old_close()
    
    # dynamically replace close function
    stream.close = new_close
'''


# TODO: add code for multiple close() calls
def open_stream(filename, mode = "r"):
    """Returns a file stream depending on the type of 'filename' and 'mode'
    
       The following types for 'filename' are handled:
       
       stream         - returns 'filename' unchanged
       iterator       - returns 'filename' unchanged
       URL string     - opens http pipe
       '-'            - opens stdin or stdout, depending on 'mode'
       other string   - opens file with name 'filename'
       
       mode is standard mode for open(): r,w,a,b
    """

    # if filename has a file interface then return it back unchanged
    if hasattr(filename, "read") or \
       hasattr(filename, "write"):
        return filename

    # if mode is reading and filename is an iterator
    if "r" in mode and hasattr(filename, "__next__"):
        return filename

    # if filename is a string then open it
    elif isinstance(filename, str):
        # open URLs
        if filename.startswith("http://"):
            import urllib.request
            return urllib.request.urlopen(filename)

        # open stdin and stdout
        elif filename == "-":
            if "w" in mode:
                return sys.stdout
            elif "r" in mode:
                return sys.stdin
            else:
                raise Exception("stream '-' can only be opened with modes r/w")

        # open regular file
        else:
            return open(filename, mode)

    # cannot handle other types for filename
    else:
        raise Exception("unknown filename type '%s'" % type(filename))
openStream = open_stream


#=============================================================================
# Delimited files
#

class DelimReader:
    """Reads delimited files"""

    def __init__(self, filename, delim=None):
        """Constructor for DelimReader
            
           arguments:
           filename  - filename or stream to read from
           delim     - delimiting character
        """

        self.infile = open_stream(filename)
        self.delim = delim

    def __iter__(self):
        return self

    def __next__(self):
        line = next(self.infile)
        fields = self.split(line)
        return fields

    def split(self, line):
        return line.rstrip().split(self.delim)


def read_delim(filename, delim=None):
    """Read an entire delimited file into memory as a 2D list"""

    return list(DelimReader(filename, delim))
readDelim = read_delim

def write_delim(filename, data, delim="\t"):
    """Write a 2D list into a file using a delimiter"""

    out = open_stream(filename, "w")
    for line in data:
        print(delim.join(map(str, line)), file=out)
writeDelim = write_delim

#=============================================================================
# printing functions
#

def default_justify(val):
    """Return the default column justification for a value.

    Numeric types (int, float) are right-justified; everything else is left.

    Args:
        val: The value whose justification is needed.

    Returns:
        "right" for int/float values, "left" otherwise.
    """
    if isinstance(val, int) or \
       isinstance(val, float):
        return "right"
    else:
        return "left"
defaultJustify = default_justify

def default_format(val):
    """Format a value for tabular display.

    Integers are formatted with comma separators via int2pretty. Percent
    values use their own __str__. Small floats use scientific notation;
    others use 4 decimal places. Everything else uses str().

    Args:
        val: The value to format.

    Returns:
        A human-readable string representation of val.
    """
    if isinstance(val, int) and \
       not isinstance(val, bool):
        return int2pretty(val)
    elif isinstance(val, Percent):
        return str(val)
    elif isinstance(val, float):
        if abs(val) < 1e-4:
            return "%.2e" % val
        else:
            return "%.4f" % val
    else:
        return str(val)
defaultFormat = default_format

def printcols(data, width=None, spacing=1, format=defaultFormat,
              justify=defaultJustify, out=sys.stdout,
              colwidth=INF, overflow="!"):
    """Prints a list or matrix in aligned columns
        
       data    - a list or matrix
       width   - maxium number of characters per line (default: 75 for lists)
       spacing - number of spaces between columns (default: 1)
       out     - stream to print to (default: sys.stdout)
    """

    if len(data) == 0:
        return

    if isinstance(data[0], list) or \
       isinstance(data[0], tuple):
        # matrix printing has default width of unlimited
        if width == None:
            width = 100000

        mat = data
    else:
        # list printing has default width 75
        if width == None:
            width = 75

        ncols = int(width / (max(map(lambda x: len(str(x)), data))+ spacing))
        mat = list2matrix(data, ncols=ncols, bycols=True)


    # turn all entries into strings
    matstr = map2(format, mat)

    # overflow
    for row in matstr:
        for j in range(len(row)):
            if len(row[j]) > colwidth:
                row[j] = row[j][:colwidth-len(overflow)] + overflow

    # ensure every row has same number of columns
    maxcols = max(map(len, matstr))
    for row in matstr:
        if len(row) < maxcols:
            row.extend([""] * (maxcols - len(row)))


    # find the maximum width char in each column
    maxwidths = map(max, map2(len, zip(* matstr)))


    # print out matrix with whitespace padding
    for i in range(len(mat)):
        fields = []
        for j in range(len(mat[i])):
            just = justify(mat[i][j])

            if just == "right":
                fields.append((" " * (maxwidths[j] - len(matstr[i][j]))) + \
                              matstr[i][j] + \
                              (" " * spacing))
            else:
                # do left by default
                fields.append(matstr[i][j] +
                              (" " * (maxwidths[j] - len(matstr[i][j]) + spacing)))
        out.write("".join(fields)[:width] + "\n")


def list2matrix(lst, nrows=None, ncols=None, bycols=True):
    """Reshape a flat list into a 2D matrix.

    Args:
        lst: The list to reshape.
        nrows: Number of rows. Inferred from ncols if not given.
        ncols: Number of columns. Inferred from nrows if not given.
            If neither is given, a roughly square shape is used.
        bycols: If True, fill the matrix column-by-column (default True).

    Returns:
        A list of lists representing the reshaped matrix.
    """

    mat = []

    if nrows == None and ncols == None:
        nrows = int(math.sqrt(len(lst)))
        ncols = int(math.ceil(len(lst) / float(nrows)))
    elif nrows == None:
        nrows = int(math.ceil(len(lst) / float(min(ncols, len(lst)))))
    else:
        ncols = int(math.ceil(len(lst) / float(min(nrows, len(lst)))))

    for i in range(nrows):
        mat.append([])
        for j in range(ncols):
            if bycols:
                k = i + j*nrows
            else:
                k = i*ncols + j
            if k < len(lst):
                mat[-1].append(lst[k])

    return mat


def printwrap(text, width=80, prefix="", out=sys.stdout):
    """Print text with line wrapping at a fixed column width.

    Args:
        text: The string to print.
        width: Maximum number of characters per line (default 80).
            If None, print the text as a single line with no wrapping.
        prefix: String prepended to each wrapped line (default "").
        out: Output stream (default sys.stdout).
    """
    if width == None:
        out.write(text)
        out.write("\n")
        return

    pos = 0
    while pos < len(text):
        out.write(prefix)
        out.write(text[pos:pos+width])
        out.write("\n")
        pos += width



def int2pretty(num):
    """Returns a pretty-printed version of an int"""

    string = str(num)
    parts = []
    l = len(string)
    for i in range(0, l, 3):
        t = l - i
        s = t - 3
        if s < 0: s = 0
        parts.append(string[s:t])
    parts.reverse()
    return ",".join(parts)


def pretty2int(string):
    """Parses a pretty-printed version of an int into an int"""
    return int(string.replace(",", ""))


def str2bool(val):
    """Correctly converts the strings "True" and "False" to the 
       booleans True and False
    """

    if val == "True":
        return True
    elif val == "False":
        return False
    else:
        raise Exception("unknown string for bool '%s'" % val)



def print_dict(dic, key=lambda x: x, val=lambda x: x,
              num=None, cmp=cmp, order=None, reverse=False,
              spacing=4, out=sys.stdout,
              format=defaultFormat,
              justify=defaultJustify):
    """Print a dictionary as an aligned two-column table.

    Args:
        dic: Dictionary to print.
        key: Function applied to keys before printing (default identity).
        val: Function applied to values before printing (default identity).
        num: Maximum number of entries to print. Defaults to all.
        cmp: Comparison function (unused in Python 3; kept for compatibility).
        order: Key function for sorting items. If None, default sort is used.
        reverse: If True, sort in descending order (default False).
        spacing: Number of spaces between columns (default 4).
        out: Output stream (default sys.stdout).
        format: Formatting function for cell values (default default_format).
        justify: Justification function for cell values (default default_justify).
    """

    if num == None:
        num = len(dic)

    dic = mapdict(dic, key=key, val=val)
    items = list(dic.items())

    if order is not None:
        items.sort(key=order, reverse=reverse)
    else:
        items.sort(reverse=reverse)

    printcols(items[:num], spacing=spacing, out=out, format=format,
              justify=justify)
printDict = print_dict


#=============================================================================
# Parsing
#

class SafeReadIter:
    """An iterator over a file handle that stops at EOF without raising an error.

    Unlike a bare for-loop over a file, this class uses readline() and raises
    StopIteration when an empty string (EOF) is encountered.
    """
    def __init__(self, infile):
        """Initialize from an open file handle.

        Args:
            infile: An open file handle to iterate over.
        """
        self.infile = infile

    def __iter__(self):
        """Return self as the iterator."""
        return self

    def __next__(self):
        """Return the next line or raise StopIteration at EOF.

        Returns:
            Next line string from the file.

        Raises:
            StopIteration: When end of file is reached.
        """
        line = self.infile.readline()
        if line == "":
            raise StopIteration
        else:
            return line

def readWord(infile, delims = [" ", "\t", "\n"]):
    """Read the next whitespace-delimited word from a file stream.

    Args:
        infile: An open file handle to read from.
        delims: List of delimiter characters (default space, tab, newline).

    Returns:
        The next word as a string, or an empty string at EOF.
    """
    word = ""

    while True:
        char = infile.read(1)
        if char == "":
            return word
        if char not in delims:
            word += char
            break

    while True:
        char = infile.read(1)
        if char == "" or char in delims:
            return word
        word += char


def readUntil(stream, chars):
    """Read from stream until one of the given characters (or EOF) is seen.

    Args:
        stream: An open file handle.
        chars: String or iterable of stop characters.

    Returns:
        A tuple (token, char) where token is the accumulated string before
        the stop character, and char is the stop character (or "" at EOF).
    """
    token = ""
    while True:
        char = stream.read(1)
        if char in chars or char == "":
            return token, char
        token += char


def readWhile(stream, chars):
    """Read from stream while characters are in the given set.

    Args:
        stream: An open file handle.
        chars: String or iterable of accepted characters.

    Returns:
        A tuple (token, char) where token is the accumulated string of
        matching characters, and char is the first non-matching character
        (or "" at EOF).
    """
    token = ""
    while True:
        char = stream.read(1)
        if char not in chars or char == "":
            return token, char
        token += char


def skipComments(infile):
    """Yield non-comment, non-blank lines from a file.

    Args:
        infile: An iterable of lines (e.g. an open file handle).

    Yields:
        Lines that do not start with "#" and are not blank.
    """
    for line in infile:
        if line.startswith("#") or line.startswith("\n"):
            continue
        yield line


class IndentStream:
    """A write-only stream wrapper that automatically indents every line.

    Tracks a current indentation depth and prepends that many spaces to the
    start of each new line. Use indent() and dedent() to change the depth.

    Attributes:
        stream: The underlying writable stream.
        linestart: True when the next character written begins a new line.
        depth: Current indentation level in spaces.
    """

    def __init__(self, stream):
        """Initialize an IndentStream wrapping the given stream.

        Args:
            stream: A filename string or writable file object to wrap.
        """
        self.stream = open_stream(stream, "w")
        self.linestart = True
        self.depth = 0

    def indent(self, num=2):
        """Increase the indentation depth.

        Args:
            num: Number of spaces to add (default 2).
        """
        self.depth += num

    def dedent(self, num=2):
        """Decrease the indentation depth, clamped to zero.

        Args:
            num: Number of spaces to remove (default 2).
        """
        self.depth -= num
        if self.depth < 0:
            self.depth = 0

    def write(self, text):
        """Write text to the underlying stream, prepending indentation as needed.

        Args:
            text: The string to write.
        """
        lines = text.split("\n")

        for line in lines[:-1]:
            if self.linestart:
                self.stream.write(" "*self.depth)
                self.linestart = True
            self.stream.write(line + "\n")

        if len(lines) > 0:
            if text.endswith("\n"):
                self.linestart = True
            else:
                self.stream.write(" "*self.depth + lines[-1])
                self.linestart = False






#=============================================================================
# file/directory functions
#
def list_files(path, ext=""):
    """Returns a list of files in 'path' ending with 'ext'"""

    files = sorted(filter(lambda x: x.endswith(ext), os.listdir(path)))
    return [os.path.join(path, x) for x in files]
listFiles = list_files


def tempfile(path, prefix, ext):
    """Generates a a temp filename 'path/prefix_XXXXXX.ext'

    DEPRECATED: use this instead
    fd, filename = temporaryfile.mkstemp(ext, prefix)
    os.close(fd)
    """

    import tempfile
    fd, filename = tempfile.mkstemp(ext, prefix, dir=path)
    import os as _os
    _os.close(fd)

    return filename


def deldir(path):
    """Recursively remove a directory"""

    # This function is slightly more complicated because of a
    # strange behavior in AFS, that creates .__afsXXXXX files

    dirs = []

    def cleandir(arg, path, names):
        for name in names:
            filename = os.path.join(path, name)
            if os.path.isfile(filename):
                os.remove(filename)
        dirs.append(path)

    # remove files
    for dp, dn, filenames in os.walk(path): cleandir(None, dp, filenames + dn)

    # remove directories
    for i in range(len(dirs)):
        # AFS work around
        afsFiles = listFiles(dirs[-i])
        for f in afsFiles:
            os.remove(f)

        while True:
            try:
                if os.path.exists(dirs[-i]):
                    os.rmdir(dirs[-i])
            except Exception:
                continue
            break


def replace_ext(filename, oldext, newext):
    """Safely replaces a file extension new a new one"""

    if filename.endswith(oldext):
        return filename[:-len(oldext)] + newext
    else:
        raise Exception("file '%s' does not have extension '%s'" % (filename, oldext))
replaceExt = replace_ext


#=============================================================================
# sorting
#


def sortrank(lst, cmp=None, key=None, reverse=False):
    """Return the indices that would sort lst.

    Args:
        lst: The list to rank.
        cmp: Comparison function (deprecated; ignored if key is provided).
        key: A one-argument function to extract a comparison key from
            each list element (default identity).
        reverse: If True, sort in descending order (default False).

    Returns:
        A list of integer indices such that [lst[i] for i in result] is sorted.
    """
    ind = list(range(len(lst)))

    if key is None:
        ind.sort(key=lambda a: lst[a], reverse=reverse)
    else:
        ind.sort(key=lambda a: key(lst[a]), reverse=reverse)
    return ind
sortInd = sortrank


def sort_together(compare, lst, *others):
    """Sort several lists based on the sorting of 'lst'"""

    ind = sortrank(lst, compare)
    lsts = [mget(lst, ind)]

    for other in others:
        lsts.append(mget(other, ind))

    return lsts
sortTogether = sort_together

def invperm(perm):
    """Return the inverse of a permutation.

    Args:
        perm: A list of unique integers 0..n-1 representing a permutation.

    Returns:
        A list inv such that inv[perm[i]] == i for all i.
    """
    inv = [0] * len(perm)
    for i in range(len(perm)):
        inv[perm[i]] = i
    return inv
invPerm = invperm



#=============================================================================
# histograms, distributions
#

def oneNorm(vals):
    """Normalize a list of values so that they sum to 1.

    Args:
        vals: A list or iterable of numeric values.

    Returns:
        A list of values each divided by the total sum.
    """
    s = float(sum(vals))
    return [x/s for x in vals]


def bucketSize(array, ndivs=None, low=None, width=None):
    """Determine bucket parameters for dividing array values into bins.

    Exactly one of ndivs or width should be supplied (or neither, which
    defaults to ndivs=20). The other value is derived from the data.

    Args:
        array: A sequence of numeric values.
        ndivs: Desired number of bins. Derived from width if not given.
        low: Lower bound for binning. Defaults to min(array).
        width: Desired bin width. Derived from ndivs if not given.

    Returns:
        A tuple (ndivs, low, width) with all three values resolved.
    """

    if low is None:
        low = min(array)

    if ndivs is None:
        if width is None:
            ndivs = 20
        else:
            ndivs = int(math.ceil(max((max(array) - low) / float(width), 1)))

    if width is None:
        width = (max(array) - low) / float(ndivs)

    return ndivs, low, width


def bucketBin(item, ndivs, low, width):
    """
    Return the bin for an item
    """

    assert item >= low, Exception("negative bucket index")
    return min(int((item - low) / width), ndivs-1)


def bucket(array, ndivs=None, low=None, width=None, key=lambda x: x):
    """Group elements of array into ndivs buckets.

    Args:
        array: A sequence of items to bucket.
        ndivs: Number of buckets (inferred if not given).
        low: Lower bound for the first bucket (default min of key values).
        width: Bucket width (inferred if not given).
        key: Function to extract a numeric comparison key from each item
            (default identity).

    Returns:
        A tuple (x, h) where x is a list of bucket lower-bound values and
        h is a list of lists containing the array items in each bucket.
    """

    keys = map(key, array)

    # set bucket sizes
    ndivs, low, width = bucketSize(keys, ndivs, low, width)

    # init histogram
    h = [[] for i in range(ndivs)]
    x = []

    # bin items
    for i in array:
        if i >= low:
            h[bucketBin(key(i), ndivs, low, width)].append(i)
    for i in range(ndivs):
        x.append(i * width + low)
    return (x, h)


def hist(array, ndivs=None, low=None, width=None):
    """Create a histogram of array values.

    Args:
        array: A sequence of numeric values.
        ndivs: Number of histogram bins (default 20 if width is also None).
        low: Lower bound of the first bin. Defaults to min(array).
        width: Bin width (inferred from ndivs if not given).

    Returns:
        A tuple (x, h) where x is a list of bin lower-bound values and
        h is a list of integer counts for each bin.
    """

    # set bucket sizes
    ndivs, low, width = bucketSize(array, ndivs, low, width)

    # init histogram
    h = [0] * ndivs
    x = []

    # count items
    for i in array:
        if i >= low:
            j = bucketBin(i, ndivs, low, width)
            if j < ndivs:
                h[j] += 1
    for i in range(ndivs):
        x.append(i * width + low)
    return (x, h)


def hist2(array1, array2,
          ndivs1=None, ndivs2=None,
          low1=None, low2=None,
          width1=None, width2=None):
    """Perform a 2D histogram over two arrays.

    Args:
        array1: First sequence of numeric values (mapped to columns).
        array2: Second sequence of numeric values (mapped to rows).
        ndivs1: Number of bins for array1 (default derived from data).
        ndivs2: Number of bins for array2 (default derived from data).
        low1: Lower bound for array1 bins. Defaults to min(array1).
        low2: Lower bound for array2 bins. Defaults to min(array2).
        width1: Bin width for array1 (inferred if not given).
        width2: Bin width for array2 (inferred if not given).

    Returns:
        A tuple (labels, h) where labels is a 2D list of [x, y] bin
        coordinates and h is a 2D list of integer counts.
    """


    # set bucket sizes
    ndivs1, low1, width1 = bucketSize(array1, ndivs1, low1, width1)
    ndivs2, low2, width2 = bucketSize(array2, ndivs2, low2, width2)

    # init histogram
    h = [[0] * ndivs1 for i in range(ndivs2)]
    labels = []

    for j,i in zip(array1, array2):
        if j > low1 and i > low2:
            h[bucketBin(i, ndivs2, low2, width2)] \
             [bucketBin(j, ndivs1, low1, width1)] += 1

    for i in range(ndivs2):
        labels.append([])
        for j in range(ndivs1):
            labels[-1].append([j * width1 + low1,
                               i * width2 + low2])
    return labels, h


def histbins(bins):
    """Convert bin start positions to bin center positions for GNUPLOT plotting.

    Args:
        bins: A list of bin start positions.

    Returns:
        A list of bin center positions the same length as bins.
    """

    bins2 = []

    if len(bins) == 1:
        bins2 = [bins[0]]
    else:
        for i in range(len(bins) - 1):
            bins2.append((bins[i] + bins[i+1]) / 2.0)
        bins2.append(bins[-1] + (bins[-1] - bins[-2]) / 2.0)

    return bins2


def distrib(array, ndivs=None, low=None, width=None):
    """Compute the probability density distribution of array.

    Normalises histogram counts by the total number of items and bin width,
    giving an approximate PDF.

    Args:
        array: A sequence of numeric values.
        ndivs: Number of bins (default derived from data).
        low: Lower bound of the first bin. Defaults to min(array).
        width: Bin width (inferred if not given).

    Returns:
        A tuple (x, h) where x is bin lower-bound values and h is a list
        of density values (count / total / width).
    """

    # set bucket sizes
    ndivs, low, width = bucketSize(array, ndivs, low, width)

    h = hist(array, ndivs, low, width)

    total = float(sum(h[1]))
    return (h[0], [(x/total)/width for x in h[1]])


def hist_int(array):
    """Returns a histogram of integers as a list of counts"""

    hist = [0]  * (max(array) + 1)
    negative = []
    for i in array:
        if (i >= 0):
            hist[i] += 1
        else:
            negative.append(i)
    return hist
histInt = hist_int

def hist_dict(array):
    """Returns a histogram of any items as a dict.
    
       The keys of the returned dict are elements of 'array' and the values
       are the counts of each element in 'array'.
    """

    hist = {}
    for i in array:
        if i in hist:
            hist[i] += 1
        else:
            hist[i] = 1
    return hist
histDict = hist_dict


def print_hist(array, ndivs=20, low=None, width=None,
              cols=75, spacing=2, out=sys.stdout):
    """Print a text-based histogram with ASCII bar chart.

    Args:
        array: A sequence of numeric values to histogram.
        ndivs: Number of bins (default 20).
        low: Lower bound for the first bin. Defaults to min(array).
        width: Bin width (inferred if not given).
        cols: Total character width of the output including bars (default 75).
        spacing: Number of spaces between columns (default 2).
        out: Output stream (default sys.stdout).
    """
    data = list(hist(array, ndivs, low=low, width=width))

    # find max bar
    maxwidths = map(max, map2(compose(len, str), data))
    maxbar = cols- sum(maxwidths) - 2 * spacing

    # make bars
    bars = []
    maxcount = max(data[1])
    for count in data[1]:
        bars.append("*" * int(count * maxbar / float(maxcount)))
    data.append(bars)

    printcols(zip(* data), spacing=spacing, out=out)
printHist = print_hist




# NOTE: rasmus library imports removed — rasmus is not Python 3 compatible.




