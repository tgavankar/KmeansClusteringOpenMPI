import sys
import csv
import random
import numpy
import getopt
import math


def usage():
    print '$> python dnaStrands.py <required args> [optional args]\n' + \
        '\t-c <#>\t\tNumber of clusters to generate\n' + \
        '\t-p <#>\t\tNumber of points per cluster (0 for random 1-100)\n' + \
        '\t-o <file>\tFilename for the output of the raw data\n' + \
        '\t-v [#]\t\tLength of DNA strand (defaults to 10)\n'


def stringDistance(s1, s2):
    count = 0
    for i in range(0, len(s1)):
        if s1[i] is not s2[i]:
            count += 1
    return count


def tooClose(string, strings, minDist):
    '''
    Computes the string distance between the string and all stringss
    in the list, and if any stringss in the list are closer than minDist,
    this method returns true.
    '''
    for pair in strings:
        if stringDistance(string, pair) < minDist:
                return True

    return False


def handleArgs(args):
    # set up return values
    numClusters = -1
    numPoints = -1
    output = None
    maxValue = 10

    try:
        optlist, args = getopt.getopt(args[1:], 'c:p:v:o:')
    except getopt.GetoptError, err:
        print str(err)
        usage()
        sys.exit(2)

    for key, val in optlist:
        # first, the required arguments
        if   key == '-c':
            numClusters = int(val)
        elif key == '-p':
            numPoints = int(val)
            if numPoints == 0:
                numPoints = random.rand(1, 100)
        elif key == '-o':
            output = val
        # now, the optional argument
        elif key == '-v':
            maxValue = int(val)

    # check required arguments were inputted
    if numClusters < 0 or numPoints < 0 or \
            maxValue < 1 or \
            output is None:
        usage()
        sys.exit()
    return (numClusters, numPoints, output, \
            maxValue)


def drawOrigin(maxLength):
    return (''.join(random.choice(['a', 't', 'g', 'c'])
            for i in xrange(maxLength)))


def mutateString(orig, distance):
    origlist = list(orig)
    toChange = random.sample([i for i in range(0, len(origlist))], distance)
    for i in toChange:
        dna = ['a', 't', 'g', 'c']
        dna.remove(origlist[i])
        origlist[i] = random.choice(dna)
    return ''.join(origlist)


def main():
    # start by reading the command line
    numClusters, \
    numPoints, \
    output, \
    maxValue = handleArgs(sys.argv)

    writer = csv.writer(open(output, "w"))

    # step 1: generate each DNA centroid
    centroids_radii = []
    minDistance = 0
    for i in range(0, numClusters):
        centroid_radius = drawOrigin(maxValue)
        # is it far enough from the others?
        while (tooClose(centroid_radius, centroids_radii, minDistance)):
            centroid_radius = drawOrigin(maxValue)
        centroids_radii.append(centroid_radius)

    # step 2: generate the strands for each centroid
    points = []
    minClusterVar = 1
    maxClusterVar = maxValue * 0.25
    for i in range(0, numClusters):
        # compute the variance for this cluster
        variance = numpy.random.uniform(minClusterVar, maxClusterVar)
        cluster = centroids_radii[i]
        for j in range(0, numPoints):
            # generate a Hamming Distance with specified variance
            # Distance is normally-distributed around centroids[i]
            distance = int(abs(numpy.random.normal(0, variance)))
            if distance > maxValue:
                distance = maxValue
            # write the strings out by mutating them based on the distance.
            writer.writerow([mutateString(cluster, distance)])


if __name__ == "__main__":
    main()
