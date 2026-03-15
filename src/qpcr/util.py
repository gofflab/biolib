'''
Miscellaneous utility functions for the qpcr package.

@author: lgoff
'''

#Misc Tools and Utilities
def uniqify(seq):
    """Return a list of unique elements from a sequence.

    Deduplicates by inserting elements into a dict. The returned order is
    not guaranteed to be the same as the input order (depends on dict
    insertion-order behaviour of the Python version).

    Args:
        seq: Any iterable of hashable elements.

    Returns:
        A list containing each unique element from ``seq`` exactly once.
    """
    # Not order preserving
    keys = {}
    for e in seq:
        keys[e] = 1
    return list(keys.keys())
