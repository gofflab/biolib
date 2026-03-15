'''
Custom data type implementations for seqlib/RNASeq data processing.

Provides a Stack (LIFO), a binary search tree (BinaryTree / BinaryNode /
EmptyNode), and a directed-graph (Graph) useful for path-finding in
acyclic graphs.

Created on Dec 14, 2009

My custom data types to help with RNASeq data

@author: lgoff
'''
error = 'stack2.error'

class Stack:
    '''
    A last-in, first-out (LIFO) stack data structure backed by a Python list.

    Supports push, pop, and peek operations, and delegates unknown attribute
    lookups to the underlying list so list methods are accessible directly.
    '''

    def __init__(self, start=[]):
        '''
        Initialize the Stack, optionally pre-loading it with items.

        Items from start are pushed in order and then reversed so that the
        first element of start ends up at the top of the stack.

        Args:
            start: An optional list of items to pre-load. Defaults to [].
        '''
        self.stack = []
        for x in start: self.push(x)
        self.reverse()

    def push(self, obj):
        """Push an item onto the top of the stack.

        Args:
            obj: The object to place on top of the stack.
        """
        self.stack = [obj] + self.stack

    def pop(self):
        """Remove and return the item at the top of the stack.

        Returns:
            The top item of the stack.

        Raises:
            stack2.error: If the stack is empty (underflow).
        """
        if not self.stack: raise error('underflow')
        top, self.stack = self.stack[0], self.stack[1:]
        return top

    def top(self):
        """Return the top item without removing it.

        Returns:
            The item currently at the top of the stack.

        Raises:
            stack2.error: If the stack is empty (underflow).
        """
        if not self.stack: raise error('underflow')
        return self.stack[0]

    def empty(self):
        """Return True if the stack contains no items.

        Returns:
            True if the stack is empty, False otherwise.
        """
        return not self.stack

    #Overloads
    def __repr__(self):
        """Return a string representation of the stack."""
        return '[Stack:%s]' % self.stack

    def __cmp__(self, other):
        """Compare this stack to another by their underlying lists."""
        return cmp(self.stack, other.stack)

    def __len__(self):
        """Return the number of items in the stack."""
        return len(self.stack)

    def __add__(self, other):
        """Concatenate two stacks and return a new Stack.

        Args:
            other: Another Stack instance to append.

        Returns:
            A new Stack containing items from this stack followed by other's.
        """
        return Stack(self.stack+other.stack)

    def __mul__(self, reps):
        """Repeat the stack contents reps times and return a new Stack.

        Args:
            reps: An integer number of times to repeat.

        Returns:
            A new Stack with the contents repeated reps times.
        """
        return Stack(self.stack * reps)

    def __getitem__(self, offset):
        """Return the item at the given index.

        Args:
            offset: An integer index into the underlying list.

        Returns:
            The item at position offset.
        """
        return self.stack[offset]

    def __getslice__(self, low, high):
        """Return a new Stack containing the slice from low to high.

        Args:
            low: The start index of the slice.
            high: The end index of the slice (exclusive).

        Returns:
            A new Stack containing the sliced elements.
        """
        return Stack(self.stack[low:high])

    def __getattr__(self, name):
        """Delegate attribute lookup to the underlying list.

        Args:
            name: The attribute name to look up on the underlying list.

        Returns:
            The attribute from the underlying list.
        """
        return getattr(self.stack, name)
        

##################
#Binary Trees
##################
class BinaryTree:
    """A binary search tree that delegates to recursive BinaryNode/EmptyNode objects.

    Stores values in sorted order and supports O(log n) average-case lookup
    and insertion. Duplicate values are silently ignored.
    """

    def __init__(self):
        """Initialize an empty BinaryTree."""
        self.tree = EmptyNode()

    def __repr__(self):
        """Return a parenthesized string representation of the tree."""
        return repr(self.tree)

    def lookup(self, value):
        """Return 1 if value exists in the tree, 0 otherwise.

        Args:
            value: The value to search for.

        Returns:
            1 if the value is present, 0 if it is not.
        """
        return self.tree.lookup(value)

    def insert(self, value):
        """Insert value into the tree, maintaining sort order.

        If value already exists in the tree it is not inserted again.

        Args:
            value: The value to insert.
        """
        self.tree = self.tree.insert(value)

class EmptyNode:
    """Sentinel node representing an empty position in a BinaryTree.

    Acts as the leaf terminator: lookup always fails and insert creates a
    new BinaryNode at this position.
    """

    def __repr__(self):
        """Return '*' to represent an empty node."""
        return "*"

    def lookup(self, value):
        """Return 0 because an empty node contains no value.

        Args:
            value: The value being searched for (unused).

        Returns:
            Always 0 (not found).
        """
        return 0

    def insert(self, value):
        """Create a new BinaryNode at this position with value.

        Args:
            value: The value to store in the new node.

        Returns:
            A new BinaryNode with empty left and right children.
        """
        return BinaryNode(self, value, self)  #Add new node at bottom
    
class BinaryNode:
    """An internal node of a binary search tree holding a value and two subtrees.

    Attributes:
        data: The value stored at this node.
        left: The left child node (values less than data).
        right: The right child node (values greater than data).
    """

    def __init__(self, left, value, right):
        """Initialize a BinaryNode with left subtree, a value, and right subtree.

        Args:
            left: The left child (a BinaryNode or EmptyNode).
            value: The value to store at this node.
            right: The right child (a BinaryNode or EmptyNode).
        """
        self.data, self.left, self.right = value, left, right

    def lookup(self, value):
        """Search for value in the subtree rooted at this node.

        Args:
            value: The value to search for.

        Returns:
            1 if value is found in this subtree, 0 otherwise.
        """
        if self.data == value:
            return 1
        elif self.data > value:
            return self.left.lookup(value)
        else:
            return self.right.lookup(value)

    def insert(self, value):
        """Insert value into the subtree rooted at this node.

        Traverses left if value is less than this node's data, right if
        greater. Equal values are ignored (no duplicates stored).

        Args:
            value: The value to insert.

        Returns:
            This node (possibly with an updated child subtree).
        """
        if self.data > value:
            self.left = self.left.insert(value)
        elif self.data < value:
            self.right = self.right.insert(value)
        return self

    def __repr__(self):
        """Return a parenthesized inorder string representation of this subtree."""
        return '( %s, %s, %s )' % (repr(self.left), repr(self.data), repr(self.right))

################
#Directed Acyclic Graphs
################
class Graph:
    """A node in a directed acyclic graph (DAG) that supports path-finding.

    Each Graph node has a label, optional data payload, and a list of
    outgoing edges to other Graph nodes. Multiple paths between nodes are
    found via depth-first search and stored as class-level state in
    Graph.solns.

    Attributes:
        name: A string label identifying this node.
        data: An optional data payload associated with this node.
        edges: A list of Graph nodes reachable from this node.
    """

    def __init__(self, label, extra=None):
        """Initialize a Graph node with a label and optional data.

        Args:
            label: A string name for this node.
            extra: An optional data object to associate with the node.
                Defaults to None.
        """
        self.name = label
        self.data = extra
        self.edges = []

    def __repr__(self):
        """Return the node's label as its string representation."""
        return self.name

    def search(self, goal):
        """Find all acyclic paths from this node to goal.

        Resets Graph.solns, performs a depth-first search via generate(),
        and sorts found paths by length (shortest first).

        Args:
            goal: A Graph node to search for.

        Returns:
            A list of paths (each path is a list of Graph nodes) from this
            node to goal, sorted by path length ascending.
        """
        Graph.solns = []
        self.generate([self], goal)
        Graph.solns.sort(lambda x, y: cmp(len(x), len(y)))
        return Graph.solns

    def generate(self, path, goal):
        """Recursively explore paths from this node towards goal.

        Appends the current path to Graph.solns when goal is reached.
        Avoids cycles by checking whether each neighbor is already in the
        current path before recursing.

        Args:
            path: A list of Graph nodes representing the current path from
                the search origin to this node.
            goal: A Graph node to find.
        """
        if self == goal:
            Graph.solns.append(path)
        else:
            for edge in self.edges:
                if edge not in path:
                    edge.generate(path + [edge], goal)
"""
def tests(Graph):
    for name in "ABCDEFG":
        exec "%s = Graph('%s')" % (name,name)
    
    A.edges = [B,E,G]
    B.edges = [C]
    C.edges = [D,E]
    D.edges = [F]
    E.edges = [C,F,G]
    G.edges = [A]
    
    A.search(G)
    for (start,stop) in [(E,D), (A,G), (G,F), (B,A), (D,A)]:
        print start.search(stop)
"""    