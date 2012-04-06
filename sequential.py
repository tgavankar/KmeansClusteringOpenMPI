import csv
import getopt
import math
import random
import sys


def euclideanDistance(p, q):
    return math.sqrt(sum([(q[i] - p[i]) ** 2 for i in range(len(p))]))


def hammingDistance(s1, s2):
    count = 0
    for i in range(0, len(s1)):
        if s1[i] is not s2[i]:
            count += 1
    return count


def point2str(p, c):
    return p


def dna2str(d, c):
    return (''.join(d), hammingDistance(d, c))


class Cluster:
    def __init__(self, points):
        self.points = points
        self.centroid = self.calcCentroid()

    def add(self, point):
        self.points.append(point)

    def clear(self):
        self.points = []

    def calcCentroid(self):
        # Perform point-wise averaging to calculate centroid for 2D points.
        if (isinstance(self.points[0][0], float)
            or isinstance(self.points[0][0], int)):
            return [(sum([p[i] for p in self.points], 0.0) /
                    float(len(self.points))) for i in
                    range(len(self.points[0]))]

        # Begin DNA String centroid calculation
        centroid = []
        randIndex = []
        countPicked = [0 for x in range(len(self.points))]
        # Iterate over the characters in the string.
        for i in range(len(self.points[0])):
            aCount = 0
            tCount = 0
            gCount = 0
            cCount = 0
            picked = 'r'
            # Iterate over all the strings in the cluster.
            for j in range(len(self.points)):
                # Increment the counter for a letter hit.
                if self.points[j][i] == 'a':
                    aCount += 1
                elif self.points[j][i] == 'g':
                    gCount += 1
                elif self.points[j][i] == 'c':
                    cCount += 1
                elif self.points[j][i] == 't':
                    tCount += 1
            # Mark the highest occuring letter in the group.
            if aCount > gCount and aCount > cCount and aCount > tCount:
                centroid.append('a')
                picked = 'a'
            elif gCount > aCount and gCount > tCount and gCount > cCount:
                centroid.append('g')
                picked = 'g'
            elif cCount > aCount and cCount > gCount and cCount > tCount:
                centroid.append('c')
                picked = 'c'
            elif tCount > aCount and tCount > gCount and tCount > cCount:
                centroid.append('t')
                picked = 't'
            else:
                # There was a tie in highest counts.
                centroid.append('r')
                maxVal = max([aCount, gCount, tCount, cCount])
                tied = []
                if aCount == maxVal:
                    tied.append('a')
                if cCount == maxVal:
                    tied.append('c')
                if gCount == maxVal:
                    tied.append('g')
                if tCount == maxVal:
                    tied.append('t')
                randIndex.append({'index': i, 'vals': tied})
            # Count the number of times the strings' character was chosen.
            for j in range(len(self.points)):
                if self.points[j][i] == picked:
                    countPicked[j] += 1
        # Pick the char that helps the string that "won" the least.
        for x in randIndex:
            tiedIndex = []
            for j in range(len(self.points)):
                for nuec in x['vals']:
                    if self.points[j][x['index']] == nuec:
                        tiedIndex.append(j)
            min = countPicked[tiedIndex[0]]
            minindex = 0
            for i in tiedIndex:
                if countPicked[i] < min:
                    min = countPicked[i]
                    minindex = i
            centroid[x['index']] = self.points[minindex][x['index']]

        self.centroid = centroid
        return centroid


def usage():
    print '$> python sequential.py <required args>\n' + \
        '\t-t <type>\t\tType of input data ("plot" or "dna")\n' + \
        '\t-k <#>\t\tNumber of clusters\n' + \
        '\t-u <#>\t\tNumber of K-means iterations\n' + \
        '\t-i <file>\tInput filename for the raw data\n' + \
        '\t-o <file>\tOutput filename\n'


def handleArgs(args):
    type = "plot"
    k = 2
    u = 0.0001
    input = None
    output = None

    if len(args) <= 1:
        usage()
        sys.exit(2)

    try:
        optlist, args = getopt.getopt(args[1:], 't:k:u:i:o:')
    except getopt.GetoptError, err:
        print str(err)
        usage()
        sys.exit(2)

    for key, val in optlist:
        if key == '-t':
            if val != "plot" and val != "dna":
                usage()
                sys.exit(2)
            type = val
        elif key == '-k':
            k = int(val)
        elif key == '-u':
            u = float(val)
        elif key == '-i':
            input = val
        elif key == '-o':
            output = val

    return (type, k, u, input, output)


def kmeanpp(X, k, distance):
    """Perform K-Mean++ initial centroid selection."""
    # Source: http://yongsun.me/2008/10/k-means-and-k-means-with-python/
    ntries = int(2 + math.log(k))
    n = len(X)
    centers = [X[random.randint(0, n)]]
    D = [distance(x, centers[0]) ** 2 for x in X]
    Dsum = reduce(lambda x, y: x + y, D)
    for _ in range(k - 1):
        bestDsum = bestIdx = -1
        for _ in range(ntries):
            randVal = random.random() * Dsum
            for i in range(n):
                if randVal <= D[i]:
                    break
                else:
                    randVal -= D[i]
            tmpDsum = reduce(lambda x, y: x + y,
                             (min(D[j], distance(X[j], X[i]) ** 2)
                             for j in xrange(n)))
            if bestDsum < 0 or tmpDsum < bestDsum:
                bestDsum, bestIdx = tmpDsum, i
        Dsum = bestDsum
        centers.append(X[bestIdx])
        D = [min(D[i], distance(X[i], X[bestIdx]) ** 2) for i in xrange(n)]
    return centers


def quickSelect(data, n):
    """Find the nth rank ordered element (the least value has rank 0)."""
    # Source: http://code.activestate.com/recipes/269554/
    data = list(data)
    if not 0 <= n < len(data):
        raise ValueError('not enough elements for the given rank')
    while True:
        pivot = random.choice(data)
        pcount = 0
        under, over = [], []
        uappend, oappend = under.append, over.append
        for elem in data:
            if elem < pivot:
                uappend(elem)
            elif elem > pivot:
                oappend(elem)
            else:
                pcount += 1
        if n < len(under):
            data = under
        elif n < len(under) + pcount:
            return pivot
        else:
            data = over
            n -= len(under) + pcount


def main():
    (type, k, u, infile, outfile) = handleArgs(sys.argv)

    if infile is None or outfile is None:
        usage()
        sys.exit(2)

    f = open(infile, "r")

    if type == "dna":
        distance = hammingDistance
        # Read strings as character lists into points from input file.
        points = [list(i.strip()) for i in f]
        numPoints = len(points)
        # Use QuickSelect to pick k distant starting centroids.
        clusters = [Cluster([c]) for c in
                    [quickSelect(points, i * (numPoints / k)
                     + numPoints / (2 * k))
                     for i in xrange(k)]]
    else:
        distance = euclideanDistance
        # Read file into lists for each coordinate.
        points = [[float(x) for x in i.split(",")] for i in f]
        numPoints = len(points)
        # Use K-Mean++ to pick k distant starting centroids.
        clusters = [Cluster([c]) for c in kmeanpp(points, k, distance)]

    while True:
        oldC = []
        # Save old cluster centroids for future comparison.
        for c in clusters:
            oldC.append(c.centroid)
            c.clear()

        # Find nearest centroid to each point and add it to cluster.
        for p in points:
            minDist = float("inf")
            for c in clusters:
                dist = distance(p, c.centroid)
                if dist < minDist:
                    minDist = dist
                    minClus = c
            minClus.add(p)

        # Calculate all the new centroids.
        for c in clusters:
            c.calcCentroid()

        # If the change in centroids is <= u, break.
        if max([distance(oldC[i], clusters[i].centroid)
                for i in range(0, len(clusters))]) <= u:
            break

    # Save output to CSV file.

    # Format for 2D points is 2 columns per cluster (x, y) with 2k total cols.
    # The first entry is the calculated centroid of each cluster.

    # Format for DNA points is 2 cols per cluster (DNA, hamming to centroid)
    # with a total of 2k cols.
    # The first entry is the calculated centroid of each cluster.

    csvOutput = []
    writer = csv.writer(open(outfile, "w"))

    if type == "dna":
        p2str = dna2str
    else:
        p2str = point2str

    for c in clusters:
        centroid = c.centroid
        csvOutput.append([])
        currClusterList = csvOutput[-1]
        currClusterList.append(p2str(centroid, centroid))
        currappend = currClusterList.append

        for p in c.points:
            currappend(p2str(p, centroid))

    maxClusterSize = max([len(c) for c in csvOutput])
    for p in xrange(maxClusterSize):
        e1 = [list(c[p]) if len(c) > p else [None, None] for c in csvOutput]
        writer.writerow(sum(e1, []))


if __name__ == "__main__":
    # Use cProfile to profile process.
    # Can be replaced with just a call to main() to disable profiling.
    import cProfile
    cProfile.run('main()')
