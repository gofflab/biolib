# python libs
"""Algorithmic data structures and search utilities for sequence analysis.

Provides Union-Find disjoint set, QuadTree spatial indexing, and binary search
implementations used throughout the seqlib package.
"""


#=============================================================================

class UnionFind:
    """An implementation of the UNION/FIND algorithm for disjoint sets.

    Supports efficient union and membership queries using path compression.
    Each UnionFind instance represents a single set; sets can be merged via
    union() and queried for shared membership via same().
    """

    def __init__(self, items):
        """Initialize a new UnionFind set containing the given items.

        Args:
            items: An iterable of hashable items to populate the initial set.
        """
        self.parent = None
        self.items = dict.fromkeys(items, 1)

    def __contains__(self):
        """Return True if item is a member of the root set."""
        return item in self.root().items

    def __len__(self):
        """Return the number of items in the root set."""
        return len(self.root().items)

    def __iter__(self):
        """Iterate over the items in the root set."""
        return iter(self.root().items)


    def add(self, item):
        """Add an item to the root set.

        Args:
            item: A hashable item to add to the set.
        """
        self.root().items[item] = 1

    def root(self):
        """Return the root UnionFind node for this set, applying path compression.

        Traverses parent pointers to find the canonical representative of the
        set. As a side effect, compresses the path by pointing this node
        directly at the root.

        Returns:
            The root UnionFind node representing this disjoint set.
        """
        node = self
        while node.parent:
            node = node.parent
        if node != self:
            self.parent = node
        return node

    def same(self, other):
        """Return True if this set and other share the same root (are in the same set).

        Args:
            other: Another UnionFind instance to compare against.

        Returns:
            True if both instances belong to the same disjoint set, False otherwise.
        """
        return self.root() == other.root()

    def union(self, other):
        """Merge this set with other so that all members belong to a single set.

        If both sets already share the same root, this is a no-op. Otherwise,
        all items from other's root are merged into this set's root, and
        other's root is reparented.

        Args:
            other: Another UnionFind instance to merge with this set.
        """
        root1 = self.root()
        root2 = other.root()
        if root1 == root2:
            return

        root1.items.update(root2.items)
        root2.items = {}
        root2.parent = root1

    def members(self):
        """Return a view of all items belonging to this set.

        Returns:
            A dict_keys view of all items in the root set.
        """
        return self.root().items.keys()


    # old function DON'T USE

    def has(self, item):
        """DEPRECATED: use x in set"""
        return item in self.members()

    def size(self):
        """DEPRECATED: use len(set)"""
        return len(self.root().items)


#=============================================================================
# QuadTree data structure

class Rect:
    """A representation of an axis-aligned rectangle.

    Stores the bounding box as (x1, y1) lower-left and (x2, y2) upper-right
    corners, normalizing the coordinates so that x1 <= x2 and y1 <= y2
    regardless of the order the arguments are supplied.
    """

    def __init__(self, x1, y1, x2, y2):
        """Initialize a Rect, normalizing so that (x1, y1) is the lower-left corner.

        Args:
            x1: X coordinate of one horizontal boundary.
            y1: Y coordinate of one vertical boundary.
            x2: X coordinate of the other horizontal boundary.
            y2: Y coordinate of the other vertical boundary.
        """
        if x1 < x2:
            self.x1 = x1
            self.x2 = x2
        else:
            self.x1 = x2
            self.x2 = x1
        if y1 < y2:
            self.y1 = y1
            self.y2 = y2
        else:
            self.y1 = y2
            self.y2 = y1

class QuadNode:
    """A single entry stored in a QuadTree leaf node.

    Associates an arbitrary item with the bounding Rect used for spatial
    indexing inside the QuadTree.
    """

    item = None
    rect = None

    def __init__(self, item, rect):
        """Initialize a QuadNode with an item and its bounding rectangle.

        Args:
            item: The object to store (any type).
            rect: A Rect instance representing the spatial extent of item.
        """
        self.item = item
        self.rect = rect


class QuadTree:
    """A spatial index that partitions 2-D space into four quadrants recursively.

    Items are stored alongside their bounding Rect. When a leaf node exceeds
    MAX items and has not yet reached MAX_DEPTH, it is split into four child
    QuadTree nodes and its items are redistributed. Items whose bounding
    rectangles span multiple quadrants are stored in every overlapping child.

    Class attributes:
        MAX: Maximum number of items in a leaf before splitting (default 10).
        MAX_DEPTH: Maximum recursion depth allowed for splits (default 10).
    """

    MAX = 10
    MAX_DEPTH = 10

    def __init__(self, x, y, size, depth = 0):
        """Initialize a QuadTree node centered at (x, y) with a given half-size.

        Args:
            x: X coordinate of this node's center.
            y: Y coordinate of this node's center.
            size: Half-width (and half-height) of the region covered by this node.
            depth: Current depth of this node in the tree (0 for the root).
        """
        self.nodes = []
        self.children = []
        self.center = [x, y]
        self.size = size
        self.depth = depth

    def insert(self, item, rect):
        """Insert an item with the given bounding rectangle into the tree.

        If this node is a leaf, the item is appended to the local node list.
        If the leaf then exceeds MAX items and depth allows, the node is split.
        If this node already has children, the item is forwarded to the
        appropriate child(ren).

        Args:
            item: The object to store.
            rect: A Rect instance representing the spatial extent of item.
        """
        if len(self.children) == 0:
            self.nodes.append(QuadNode(item, rect))

            if len(self.nodes) > self.MAX and self.depth < self.MAX_DEPTH:
                self.split()
        else:
            self.insertIntoChildren(item, rect)

    def insertIntoChildren(self, item, rect):
        """Forward an item into every child quadrant that its bounding rect overlaps.

        The four children are ordered: [bottom-left, top-left, bottom-right,
        top-right] relative to the center of this node.

        Args:
            item: The object to store.
            rect: A Rect instance representing the spatial extent of item.
        """
        if rect.x1 < self.center[0]:
            if rect.y1 < self.center[1]:
                self.children[0].insert(item, rect)
            if rect.y2 > self.center[1]:
                self.children[1].insert(item, rect)
        if rect.x2 > self.center[0]:
            if rect.y1 < self.center[1]:
                self.children[2].insert(item, rect)
            if rect.y2 > self.center[1]:
                self.children[3].insert(item, rect)

    def split(self):
        """Split this leaf node into four child QuadTree nodes.

        Creates four children covering the four quadrants of this node's
        region, then redistributes all currently held items into the children.
        After splitting, the local node list is cleared.
        """
        self.children = [QuadTree(self.center[0] - self.size/2,
                                  self.center[1] - self.size/2,
                                  self.size/2, self.depth + 1),
                         QuadTree(self.center[0] - self.size/2,
                                  self.center[1] + self.size/2,
                                  self.size/2, self.depth + 1),
                         QuadTree(self.center[0] + self.size/2,
                                  self.center[1] - self.size/2,
                                  self.size/2, self.depth + 1),
                         QuadTree(self.center[0] + self.size/2,
                                  self.center[1] + self.size/2,
                                  self.size/2, self.depth + 1)]

        for node in self.nodes:
            self.insertIntoChildren(node.item, node.rect)
        self.nodes = []

    def query(self, rect, results = {}, ret = True):
        """Return all items whose bounding rectangles overlap the query rect.

        Recursively traverses child nodes that overlap rect. At leaf nodes,
        items whose stored Rect intersects rect are added to the result set.
        The results dict is used for deduplication (items are keys).

        Args:
            rect: A Rect instance defining the query region.
            results: A dict used internally to accumulate results across
                recursive calls. Callers should not pass this argument.
            ret: If True (the default for the top-level call), the method
                returns the keys of the results dict. Recursive calls pass
                False to suppress the return.

        Returns:
            A dict_keys view of all items that overlap rect, or None when
            called recursively (ret=False).
        """
        if ret:
            results = {}

        if len(self.children) > 0:
            if rect.x1 < self.center[0]:
                if rect.y1 < self.center[1]:
                    self.children[0].query(rect, results, False)
                if rect.y2 > self.center[1]:
                    self.children[1].query(rect, results, False)
            if rect.x2 > self.center[0]:
                if rect.y1 < self.center[1]:
                    self.children[2].query(rect, results, False)
                if rect.y2 > self.center[1]:
                    self.children[3].query(rect, results, False)
        else:
            for node in self.nodes:
                if node.rect.x2 > rect.x1 and node.rect.x1 < rect.x2 and \
                   node.rect.y2 > rect.y1 and node.rect.y1 < rect.y2:
                    results[node.item] = True

        if ret:
            return results.keys()

    def getSize(self):
        """Return the total number of items stored in this node and all descendants.

        Returns:
            An integer count of all QuadNode items held in the subtree rooted
            at this node.
        """
        size = 0
        for child in self.children:
            size += child.getSize()
        size += len(self.nodes)
        return size

#=============================================================================
# TODO: make a funtion based linear search

def binsearch(lst, val, compare=None, order=1):
    """Performs binary search for val in lst using compare

       if val in lst:
          Returns (i, i) where lst[i] == val
       if val not in lst
          Returns index i,j where
            lst[i] < val < lst[j]

       runs in O(log n)
    """
    if compare is None:
        def compare(a, b): return (a > b) - (a < b)

    assert order == 1 or order == -1

    low = 0
    top = len(lst) - 1

    if len(lst) == 0:
        return None, None

    if compare(lst[-1], val) * order == -1:
        return (top, None)

    if compare(lst[0], val) * order == 1:
        return (None, low)

    while top - low > 1:
        ptr = (top + low) // 2

        comp = compare(lst[ptr], val) * order

        if comp == 0:
            # have we found val exactly?
            return ptr, ptr
        elif comp == -1:
            # is val above ptr?
            low = ptr
        else:
            top = ptr


    # check top and low for exact hits
    if compare(lst[low], val) == 0:
        return low, low
    elif compare(lst[top], val) == 0:
        return top, top
    else:
        return low, top



if __name__ == "__main__":

    if True:
        set1 = UnionFind()
        set2 = UnionFind()
        set3 = UnionFind()

        set1.add(1)
        set1.add(2)
        print(set1.size())
        set2.add(3)
        set2.add(4)
        set2.add(5)
        print(set2.size())
        set3.add(5)
        set3.add(6)
        set3.add(7)
        print(set3.size())
        print(set1.same(set2))
        set1.union(set2)
        print(set1.same(set2))
        set1.union(set3)

        print(set1.members())
        print(set1.size(), set2.size())
