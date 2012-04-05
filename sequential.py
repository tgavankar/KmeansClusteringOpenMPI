import csv
import getopt
import math
import random
import sys

def euclideanDistance(p, q):
    return math.sqrt(sum([(q[i]-p[i])**2 for i in range(len(p))]))

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
        if isinstance(self.points[0][0], float) or isinstance(self.points[0][0], int):
            return [sum([p[i] for p in self.points], 0.0) / float(len(self.points)) for i in range(len(self.points[0]))]
        # Begin DNA String centroid calculation
        centroid = []
        randIndex = []
        countPicked = [0 for x in range(len(self.points))]
        for i in range(len(self.points[0])):
            aCount = 0
            tCount = 0
            gCount = 0
            cCount = 0
            picked = 'r'
            for j in range(len(self.points)):
                if self.points[j][i] == 'a':
                    aCount+=1;
                elif self.points[j][i] == 'g':
                    gCount+=1;
                elif self.points[j][i] == 'c':
                    cCount+=1;
                elif self.points[j][i] == 't':
                    tCount+=1;
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
                centroid.append('r')
                maxVal = max([aCount,gCount,tCount,cCount])
                tied = []
                if aCount == maxVal:
                    tied.append('a')
                if cCount == maxVal:
                    tied.append('c')
                if gCount == maxVal:
                    tied.append('g')
                if tCount == maxVal:
                    tied.append('t')
                randIndex.append({'index':i,'vals':tied})
            for j in range(len(self.points)):
                if self.points[j][i] == picked:
                    countPicked[j]+=1;
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
        points = [list(i.strip()) for i in f]
    else:
        distance = euclideanDistance
        points = [[float(x) for x in i.split(",")] for i in f]
    
    numPoints = len(points)

    clusters = [Cluster([c]) for c in [quickSelect(points, i*(numPoints/k)+numPoints/(2*k)) for i in xrange (k)]] 

    while True:
        oldC = []
        for c in clusters:
            oldC.append(c.centroid)
            c.clear()    

        for p in points:
            minDist = float("inf")
            for c in clusters:        
                dist = distance(p, c.centroid)
                if dist < minDist:
                    minDist = dist
                    minClus = c
            minClus.add(p)

        for c in clusters:
            c.calcCentroid()
    
        if max([distance(oldC[i], clusters[i].centroid) for i in range(0, len(clusters))]) < u:
            break

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
                   
    csvTranspose = zip(*csvOutput)

    for e in csvTranspose:
        e1 = [list(i) for i in e]
        writer.writerow(sum(e1, []))


if __name__ == "__main__": 
    main()
