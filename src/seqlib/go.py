"""Gene Ontology (GO) database parsing and traversal utilities.

Provides classes and functions for loading a Gene Ontology OBO-XML file,
representing GO terms, and traversing the GO DAG to retrieve all ancestor
terms for a given GO accession.  Includes deprecated tab-delimited annotation
file readers.
"""
import xml.sax.handler
from xml.sax import make_parser
from xml.sax.handler import feature_namespaces


def readGo(filename):
    """Reads a tab-delimited GO annotation file and returns a mapping of gene IDs to GO terms.

    DEPRECATED: This function relies on the Python 2 built-in file() and the
    non-standard Dict class.  It is retained for historical reference only.

    Args:
        filename: Path to a tab-delimited GO annotation file where column 0
            contains the gene/feature identifier and column 4 contains the
            GO term.  Lines containing 'GI:' are skipped.

    Returns:
        A Dict (default list) mapping gene identifiers to lists of GO term
        strings.
    """
    terms = Dict(default=[])
    
    for line in file(filename):
        if "GI:" in line:# or "KEGG:" in line:
            continue
        tokens = line.rstrip().split("\t")
        try:
            terms[tokens[0]].append(tokens[4])
        except:
            print(line)
    
    return terms


def readCommonNames(filename):
    """Reads a tab-delimited file mapping identifiers to common gene names.

    DEPRECATED: Relies on the Python 2 built-in file().  Retained for
    historical reference only.

    Args:
        filename: Path to a two-column tab-delimited file where column 0 is
            the primary identifier and column 1 is the common name ('-'
            entries are skipped).

    Returns:
        A dictionary mapping primary identifiers to common name strings.
    """
    commonNames = {}

    for line in file(filename):
        tokens = line.rstrip().split("\t")

        if tokens[1] != '-':
            commonNames[tokens[0]] = tokens[1]
    return commonNames



class GoTerm:
    """Represents a single Gene Ontology term.

    Attributes:
        accession: GO accession string (e.g. 'GO:0008150').
        name: Human-readable term name (e.g. 'biological process').
        definition: Textual definition of the term.
        is_a: List of parent GO accession strings linked by 'is_a' relations.
        part_of: List of parent GO accession strings linked by 'part_of'
            relations.
    """
    def __init__(self):
        """Initialises a GoTerm with empty/default attribute values."""
        self.accession = ""
        self.name = ""
        self.definition = ""
        self.is_a = []
        self.part_of = []
#        self.synonym = []

class AllTerm(GoTerm):
    """Synthetic top-level GO term used as the root of the GO hierarchy.

    AllTerm has a fixed accession and name of 'all' and is added to the
    GoDatabase after parsing to provide a single root node for traversal.
    """
    def __init__(self):
        """Initialises AllTerm with accession='all' and name='all'."""
        GoTerm.__init__(self)

        self.accession = "all"
        self.name = "all"
        self.defintion = "top-level term"

class GoHandler(xml.sax.handler.ContentHandler):
    """SAX content handler for parsing Gene Ontology OBO-XML files.

    Builds a dictionary of GoTerm objects from a GO OBO-XML file as it is
    streamed through a SAX parser.  Handles go:term, go:is_a, go:part_of,
    go:accession, go:name, and go:definition elements.

    Attributes:
        terms: Dictionary mapping GO accession strings to GoTerm objects.
        term: The GoTerm currently being parsed, or None between terms.
        elm: Name of the XML element currently open, used to route character
            data to the correct GoTerm attribute.
        base: URL prefix for the GO namespace, used to strip absolute URIs
            to relative accession strings in is_a and part_of relations.
    """
    def __init__(self, base):
        """Initialises the GoHandler with a namespace base URL.

        Args:
            base: URL prefix for the GO namespace
                (e.g. 'http://www.geneontology.org/go#').
        """
        self.terms = {}
        self.term = None
        self.elm = ""
        self.base = base

    def startElement(self, name, attrs):
        """Handles the opening of an XML element during SAX parsing.

        Creates a new GoTerm when a go:term element opens, and appends
        parent accessions when go:is_a or go:part_of elements are encountered.

        Args:
            name: Local name of the XML element.
            attrs: AttributesImpl object providing element attributes.
        """
        if name == "go:term":
            self.term = GoTerm()
        elif name == "go:is_a":
            ref = attrs["rdf:resource"]
            if ref.startswith(self.base):
                self.term.is_a.append(ref[len(self.base):])
        elif name == "go:part_of":
            ref = attrs["rdf:resource"]
            if ref.startswith(self.base):
                self.term.part_of.append(ref[len(self.base):])
        self.elm = name
    
    def endElement(self, name):
        """Handles the closing of an XML element during SAX parsing.

        Stores the completed GoTerm in the terms dictionary when a go:term
        element closes, and resets the current element tracker.

        Args:
            name: Local name of the closing XML element.
        """
        if name == "go:term":
            self.terms[self.term.accession] = self.term
        self.elm = ""
    
    def characters(self, text):
        if self.elm == "go:accession":
            self.term.accession = text
        elif self.elm == "go:name":
            self.term.name = text
        elif self.elm == "go:definition":
            self.term.definition = text
        

class GoDatabase:
    def __init__(self, filename):
        # Create a parser
        parser = make_parser()

        # Tell the parser we are not interested in XML namespaces
        parser.setFeature(feature_namespaces, 0)

        # Create the handler
        dh = GoHandler("http://www.geneontology.org/go#")

        # Tell the parser to use our handler
        parser.setContentHandler(dh)

        # Parse the input
        parser.parse(filename)

        self.terms = dh.terms
        
        # add top level term
        self.terms["all"] = AllTerm()
    
    
    def getAllParents(self, goid, touched=None, count=0, ret=True):
        if touched == None:
            touched = {}
        
        if goid in self.terms:
            term = self.terms[goid]
            parents =  term.is_a + term.part_of
            
            for parent in parents:
                if parent not in touched and parent != "all":
                    touched[parent] = count
                    count += 1
            
            for parent in parents:
                self.getAllParents(parent, touched, count, False)
        
        if ret:
            parents = touched.keys()
            parents.sort(key=lambda x: touched[x])
            return parents
