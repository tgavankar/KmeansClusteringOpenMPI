import getopt
from heapq import nlargest
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
        '\t-i <file>\tFilename for the raw data\n'

def handleArgs(args):
    type = "plot"
    k = 2
    u = 0.0001
    input = None

    try:
        optlist, args = getopt.getopt(args[1:], 't:k:u:i:')
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

    return (type, k, u, input)

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

def main():
    (type, k, u, fp) = handleArgs(sys.argv)
    f = open(fp, "r")

    if type == "dna":
        distance = hammingDistance
        points = [list(i.strip()) for i in f]
    else:
        distance = euclideanDistance
        points = [[float(x) for x in i.split(",")] for i in f]
    
    numPoints = len(points)
    clusters = [Cluster([c]) for c in select(points, [i*(numPoints/k)+numPoints/(2*k) for i in range(0, k)])]
    
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

        if max([distance(oldC[i], clusters[i].centroid) for i in range(0, len(clusters))]) < u:
            break


    for c in clusters:
        print "Centroid: %s" % c.centroid
        print "Points:"
        for p in c.points:
            print "%s: %s" % (''.join(p), hammingDistance(p, c.centroid))
        

if __name__ == "__main__": 
    main()
