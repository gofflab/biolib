'''
Created on Dec 14, 2009

My custom data types to help with RNASeq data

@author: lgoff
'''
error = 'stack2.error'

class Stack:
    '''
    Basic 'stack' data type
    '''
    def __init__(self, start=[]):
        '''
        Constructor
        '''
        self.stack = []
        for x in start: self.push(x)
        self.reverse()
    
    def push(self,obj):
        self.stack = [obj] + self.stack
    
    def pop(self):
        if not self.stack: raise error, 'underflow'
        top, self.stack = self.stack[0], self.stack[1:]
        return top
    
    def top(self):
        if not self.stack: raise error, 'underflow'
        return self.stack[0]
    
    def empty(self):
        return not self.stack
    
    #Overloads
    def __repr__(self):
        return '[Stack:%s]' % self.stack
    
    def __cmp__(self,other):
        return cmp(self.stack, other.stack)
    
    def __len__(self):
        return len(self.stack)
    
    def __add__(self,other):
        return Stack(self.stack+other.stack)
    
    def __mul__(self,reps):
        return Stack(self.stack * reps)
    
    def __getitem__(self,offset):
        return self.stack[offset]
    
    def __getslice__(self,low,high):
        return Stack(self.stack[low:high])
    
    def __getattr__(self,name):
        return getattr(self.stack,name)
        

##################
#Binary Trees
##################
class BinaryTree:
    def __init__(self):
        self.tree = EmptyNode()
    def __repr__(self):
        return `self.tree`
    def lookup(self,value):
        return self.tree.lookup(value)
    def insert(self,value):
        self.tree = self.tree.insert(value)

class EmptyNode:
    def __repr__(self):
        return "*"
    def lookup(self,value):                 #Fail at the bottom
        return 0
    def insert(self,value):
        return BinaryNode(self,value,self)  #Add new node at bottom
    
class BinaryNode:
    def __init__(self,left,value,right):
        self.data,self.left,self.right = value,left,right
    def lookup(self,value):
        if self.data == value:
            return 1
        elif self.data>value:
            return self.left.lookup(value)
        else:
            return self.right.lookup(value)
    def insert(self,value):
        if self.data > value:
            self.left = self.left.insert(value)
        elif self.data < value:
            self.right = self.right.insert(value)
        return self
    def __repr__(self):
        return '( %s, %s, %s )' % (`self.left`, `self.data`, `self.right`)

################
#Directed Acyclic Graphs
################
class Graph:
    def __init__(self,label,extra=None):
        self.name = label
        self.data = extra
        self.edges = []
    def __repr__(self):
        return self.name
    def search(self,goal):
        Graph.solns = []
        self.generate([self],goal)
        Graph.solns.sort(lambda x,y: cmp(len(x), len(y)))
        return Graph.solns
    def generate(self,path,goal):
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