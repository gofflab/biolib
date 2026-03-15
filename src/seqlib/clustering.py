'''
K-means clustering implementation for n-dimensional point data.

Provides Point and Cluster data structures along with a K-means clustering
algorithm and Euclidean distance metric for grouping arbitrary numeric
coordinate data.

Created on Nov 26, 2010

@author: lgoff
'''
import math
import random
import sys


#Classes
class Point:
    """A point in n-dimensional space.

    Attributes:
        coords: A list of numeric coordinates, one per dimension.
        n: The number of dimensions (length of coords).
        reference: An optional object associated with this point (e.g. an
            original data record).
    """

    # -- The Point class represents points in n-dimensional space
    # Instance variables
    # self.coords is a list of coordinates for this Point
    # self.n is the number of dimensions this Point lives in (ie, its space)
    # self.reference is an object bound to this Point
    # Initialize new Points
    def __init__(self, coords, reference=None):
        """Initialize a Point with a coordinate list and optional reference.

        Args:
            coords: A list of numeric values representing the coordinates of
                this point in n-dimensional space.
            reference: An optional object to associate with this point.
                Defaults to None.
        """
        self.coords = coords
        self.n = len(coords)
        self.reference = reference

    # Return a string representation of this Point
    def __repr__(self):
        """Return a string representation of the coordinate list."""
        return str(self.coords)

class Cluster:
    """A cluster of Points in n-dimensional space used by the K-means algorithm.

    All Points in a Cluster must share the same number of dimensions. The
    cluster maintains a centroid (the coordinate-wise mean of its points)
    which is recalculated whenever the cluster's membership changes.

    Attributes:
        points: A list of Point objects belonging to this cluster.
        n: The number of dimensions of the Points in this cluster.
        centroid: A Point representing the sample mean of all cluster points.
    """

    # -- The Cluster class represents clusters of points in n-dimensional space
    # Instance variables
    # self.points is a list of Points associated with this Cluster
    # self.n is the number of dimensions this Cluster's Points live in
    # self.centroid is the sample mean Point of this Cluster
    def __init__(self, points):
        """Initialize a Cluster from a non-empty list of same-dimensional Points.

        Args:
            points: A non-empty list of Point objects, all with the same
                number of dimensions.

        Raises:
            Exception: If points is empty ('ILLEGAL: EMPTY CLUSTER').
            Exception: If points contain mixed dimensionality
                ('ILLEGAL: MULTISPACE CLUSTER').
        """
        # We forbid empty Clusters (they don't make mathematical sense!)
        if len(points) == 0: raise Exception("ILLEGAL: EMPTY CLUSTER")
        self.points = points
        self.n = points[0].n
        # We also forbid Clusters containing Points in different spaces
        # Ie, no Clusters with 2D Points and 3D Points
        for p in points:
            if p.n != self.n: raise Exception("ILLEGAL: MULTISPACE CLUSTER")
        # Figure out what the centroid of this Cluster should be
        self.centroid = self.calculateCentroid()

    # Return a string representation of this Cluster
    def __repr__(self):
        """Return a string representation of the list of Points in this cluster."""
        return str(self.points)

    # Update function for the K-means algorithm
    # Assigns a new list of Points to this Cluster, returns centroid difference
    def update(self, points):
        """Replace this cluster's points and return how far the centroid moved.

        Used during each iteration of the K-means algorithm to reassign points
        and measure convergence.

        Args:
            points: A new list of Point objects to assign to this cluster.

        Returns:
            The Euclidean distance between the old centroid and the new
            centroid after recalculation.
        """
        old_centroid = self.centroid
        self.points = points
        self.centroid = self.calculateCentroid()
        return getDistance(old_centroid, self.centroid)

    # Calculates the centroid Point - the centroid is the sample mean Point
    # (in plain English, the average of all the Points in the Cluster)
    def calculateCentroid(self):
        """Compute and return the centroid of the current cluster points.

        The centroid is the coordinate-wise arithmetic mean of all Points
        in the cluster.

        Returns:
            A new Point whose coordinates are the mean of each dimension
            across all points in the cluster.
        """
        centroid_coords = []
        # For each coordinate:
        for i in range(self.n):
            # Take the average across all Points
            centroid_coords.append(0.0)
            for p in self.points:
                centroid_coords[i] = centroid_coords[i]+p.coords[i]
            centroid_coords[i] = centroid_coords[i]/len(self.points)
        # Return a Point object using the average coordinates
        return Point(centroid_coords)

# -- Return Clusters of Points formed by K-means clustering
def kmeans(points, k, cutoff):
    """Cluster points into k groups using the K-means algorithm.

    Randomly selects k seed points and iteratively reassigns every point to
    the nearest cluster centroid, updating centroids after each round. Stops
    when the largest centroid shift in a single iteration falls below cutoff.

    Args:
        points: A list of Point objects to cluster. All points must have the
            same dimensionality.
        k: The number of clusters to form.
        cutoff: A float convergence threshold. Iteration stops when the
            maximum centroid displacement across all clusters is less than
            this value.

    Returns:
        A list of k Cluster objects containing the final cluster assignments.
    """
    # Randomly sample k Points from the points list, build Clusters around them
    initial = random.sample(points, k)
    clusters = []
    for p in initial: clusters.append(Cluster([p]))
    # Enter the program loop
    while True:
        # Make a list for each Cluster
        lists = []
        for c in clusters: lists.append([])
        # For each Point:
        for p in points:
            # Figure out which Cluster's centroid is the nearest
            smallest_distance = getDistance(p, clusters[0].centroid)
            index = 0
            for i in range(len(clusters[1:])):
                distance = getDistance(p, clusters[i+1].centroid)
                if distance < smallest_distance:
                    smallest_distance = distance
                    index = i+1
            # Add this Point to that Cluster's corresponding list
            lists[index].append(p)
        # Update each Cluster with the corresponding list
        # Record the biggest centroid shift for any Cluster
        biggest_shift = 0.0
        for i in range(len(clusters)):
            shift = clusters[i].update(lists[i])
            biggest_shift = max(biggest_shift, shift)
        # If the biggest centroid shift is less than the cutoff, stop
        if biggest_shift < cutoff: break
    # Return the list of Clusters
    return clusters

######
#Distance metrics
######
# -- Get the Euclidean distance between two Points
def getDistance(a, b):
    """Return the Euclidean distance between two Points.

    Args:
        a: A Point object.
        b: A Point object in the same dimensional space as a.

    Returns:
        A float representing the Euclidean (straight-line) distance between
        the two points.

    Raises:
        Exception: If a and b have different numbers of dimensions
            ('ILLEGAL: NON-COMPARABLE POINTS').
    """
    # Forbid measurements between Points in different spaces
    if a.n != b.n: raise Exception("ILLEGAL: NON-COMPARABLE POINTS")
    # Euclidean distance between a and b is sqrt(sum((a[i]-b[i])^2) for all i)
    ret = 0.0
    for i in range(a.n):
        ret = ret+pow((a.coords[i]-b.coords[i]), 2)
    return math.sqrt(ret)

###########
#Utility functions
###########
# -- Create a random Point in n-dimensional space
def makeRandomPoint(n, lower, upper):
    """Create a Point with n random coordinates drawn uniformly from [lower, upper].

    Args:
        n: The number of dimensions for the new point.
        lower: The lower bound of the uniform distribution.
        upper: The upper bound of the uniform distribution.

    Returns:
        A Point object with n coordinates each sampled from
        random.uniform(lower, upper).
    """
    coords = []
    for i in range(n): coords.append(random.uniform(lower, upper))
    return Point(coords)

##############
#Main
##############
def main(args):
    """Run a demo K-means clustering on randomly generated 2-D points.

    Creates 10 random points in 2-D space within [-200, 200] and clusters
    them into 3 groups with a convergence cutoff of 0.5, then prints the
    points and resulting clusters to stdout.

    Args:
        args: Command-line argument list (not currently used).
    """
    num_points, n, k, cutoff, lower, upper = 10, 2, 3, 0.5, -200, 200
    # Create num_points random Points in n-dimensional space
    points = []
    for i in range(num_points): points.append(makeRandomPoint(n, lower, upper))
    # Cluster the points using the K-means algorithm
    clusters = kmeans(points, k, cutoff)
    # Print the results
    print("\nPOINTS:")
    for p in points: print("P:", p)
    print("\nCLUSTERS:")
    for c in clusters: print("C:", c)

if __name__=="__main__":
    main(sys.argv)
