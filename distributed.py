import csv
import getopt
import math
import random
import sys

from mpi4py import MPI

def euclideanDistance(p, q):
    return math.sqrt(sum([(q[i]-p[i])**2 for i in range(len(p))]))

def hammingDistance(s1, s2):
    count = 0
    for i in range(0, len(s1)):
        if s1[i] is not s2[i]:
            count += 1
    return count

class Cluster:
    def __init__(self, points):
        self.points = points
        self.centroid = self.calcCentroid()   
 
    def add(self, point):
        self.points.append(point)
        self.centroid = self.calcCentroid()

    def clear(self):
        self.points = []

    def calcCentroid(self):
        if isinstance(self.points[0][0], float) or isinstance(self.points[0][0], int):
            return [sum([p[i] for p in self.points], 0.0) / float(len(self.points)) for i in range(len(self.points[0]))]
        #######
        # CODE TO GET CENTROID OF DNA GOES HERE
        # points = [['a', 'g', 'a'], ['c', 'g', 'a'], ['c', 't', 'a']]
        ######
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
        return centroid

def usage():
    print '$> python sequential.py <required args>\n' + \
        '\t-t <type>\t\tType of input data ("plot" or "dna")\n' + \
        '\t-k <#>\t\tNumber of clusters\n' + \
        '\t-u <#>\t\tNumber of K-means iterations\n' + \
        '\t-i <file>\tInput filename for the raw data\n' + \
        '\t-o <file>\tOutput filename for results\n'

def handleArgs(args):
    type = "plot"
    k = 2
    u = 0.0001
    input = None
    output = None

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

def select(data, positions, start=0, end=None):
    '''For every n in *positions* find nth rank ordered element in *data*
        inplace select'''
    # Source: http://code.activestate.com/recipes/577477-select-some-nth-smallest-elements-quickselect-inpl/
    if not end: end = len(data) - 1
    if end < start:
        return []
    if end == start:
        return [data[start]]
    pivot_rand_i = random.randrange(start,end)
    pivot_rand = data[pivot_rand_i] # get random pivot
    data[end], data[pivot_rand_i] = data[pivot_rand_i], data[end]
    pivot_i = start
    for i in xrange(start, end): # partitioning about the pivot
        if data[i] < pivot_rand:
            data[pivot_i], data[i] = data[i], data[pivot_i]
            pivot_i += 1
    data[end], data[pivot_i] = data[pivot_i], data[end]
    under_positions, over_positions, mid_positions = [],[],[]
    for position in positions:
        if position == pivot_i:
            mid_positions.append(position)
        elif position < pivot_i:
            under_positions.append(position)
        else:
            over_positions.append(position)

    result = []
    if len(under_positions) > 0:
        result.extend(select(data, under_positions, start, pivot_i-1))
    if len(mid_positions) > 0:
        result.extend([data[position] for position in mid_positions])
    if len(over_positions) > 0:
        result.extend(select(data, over_positions, pivot_i+1, end))
    return result


def chunks(lst, n, withStarts=False):
    # http://stackoverflow.com/questions/2659900/python-slicing-a-list-into-n-nearly-equal-length-partitions
    division = len(lst) / float(n)
    out = [lst[int(round(division * i)): int(round(division * (i + 1)))] for i in xrange(n)]
    if withStarts:
        out = (out, [int(round(division * i)) for i in xrange(n)])
    return out

def main():
    (type, k, u, infile, outfile) = handleArgs(sys.argv)
    f = open(infile, "r")

    if type == "dna":
        distance = hammingDistance
        points = [list(i.strip()) for i in f]
    else:
        distance = euclideanDistance
        points = [[float(x) for x in i.split(",")] for i in f]
    
    numPoints = len(points)
    clusters = [Cluster([c]) for c in select(points, [i*(numPoints/k)+numPoints/(2*k) for i in range(0, k)])]
    
    comm = MPI.COMM_WORLD

    myrank = comm.Get_rank()
    nprocs = comm.Get_size()

    while True:
        oldC = []
        if myrank == 0:
            for c in clusters:
                oldC.append(c.centroid)
                c.clear()    

        clusters = comm.bcast(clusters, root=0)

        data = comm.scatter(chunks(points, nprocs), root=0)

        reply = {}

        for p in data:
            minDist = float("inf")
            for i in range(len(clusters)):
                c = clusters[i]
                dist = distance(p, c.centroid)
                if dist < minDist:
                    minDist = dist
                    minClusI = i
            reply.setdefault(minClusI, []).append(p)

        reply = comm.gather(reply, root=0)

        if myrank == 0:
            #  MAY CONTAIN DUPLICATES IN x.get(i, [])+y.get(i, [])
            toAdd = reduce(lambda x, y: dict((i, x.get(i, [])+y.get(i, [])) for i in set(x.keys()+y.keys())), reply)
        else:
            toAdd = {}

        toAdd = comm.bcast(toAdd, root=0)

        (splitClusters, indexes) = chunks(clusters, nprocs, True)

        clusData = comm.scatter(splitClusters, root=0)
        indexes = comm.bcast(indexes, root=0)

        startIndex = indexes[myrank]

        for i in range(len(clusData)):
            clusData[i].clear()
            [clusData[i].add(p) for p in toAdd.get(startIndex+i, [])]

        gatheredClus = comm.gather(clusData, root=0)

        endReached = False
        
        if myrank == 0:
            clusters = reduce(lambda x, y: x+y, gatheredClus)
            if max([distance(oldC[i], clusters[i].centroid) for i in range(0, len(clusters))]) < u:
                endReached = True

        breakOut = comm.bcast(endReached, root=0)

        if breakOut:
            break
     
    if myrank == 0:
        writer = csv.writer(open(outfile, "w"))

        for c in clusters:
            print "Centroid: %s" % c.centroid
            writer.writerow(["Centroid"] + c.centroid)
            print "Points:"
            for p in c.points:
                if type == "dna":
                    print "%s: %s" % (''.join(p), hammingDistance(p, c.centroid))
                    writer.writerow([''.join(p)])
                else:
                    print "%s\t%s" % (p[0], p[1])
                    writer.writerow(p)
         

if __name__ == "__main__": 
    main()
