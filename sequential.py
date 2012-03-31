import getopt
import math
import random
import sys

def euclideanDistance(p, q):
    return math.sqrt(sum([(q[i]-p[i])**2 for i in range(len(p))]))

def hammingDistance(p, q):
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
        return [sum([p[i] for p in self.points], 0.0) / float(len(self.points)) for i in range(len(self.points[0]))]

def usage():
    print '$> python sequential.py <required args>\n' + \
        '\t-k <#>\t\tNumber of clusters\n' + \
        '\t-u <#>\t\tNumber of K-means iterations\n' + \
        '\t-i <file>\tFilename for the raw data\n'

def handleArgs(args):
    k = 2
    u = 0.0001
    input = None

    try:
        optlist, args = getopt.getopt(args[1:], 'k:u:i:')
    except getopt.GetoptError, err:
        print str(err)
        usage()
        sys.exit(2)

    for key, val in optlist:
        if   key == '-c':
            numClusters = int(val)
        elif key == '-u':
            u = float(val)
        elif key == '-i':
            input = val

    return (k, u, input)

def main():
    (k, u, fp) = handleArgs(sys.argv)

    distance = euclideanDistance

    f = open(fp, "r")

    points = [[float(x) for x in i.split(",")] for i in f]

    # Fully random points
    #points = [[random.uniform(0, 50), random.uniform(0, 50)] for i in range(20)]

    clusters = [Cluster([c]) for c in random.sample(points, k)]

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
            print "%s\t%s" % (p[0], p[1])
        

if __name__ == "__main__": 
    main()
