'''

'''
#==========================================================================================
import sys
import operator
#==========================================================================================
def main():
    #Check argumments.
    args = parseArgs(sys.argv)
    if args == ():
        return
    inputFile,outputFile = args

    fin = openFileReading(inputFile)
    if fin == ():
        return()

    myHash = {0: 'A', 1: 'C', 2: 'G', 3: 'T', 4: '-'}
    fout = openFileWriting(outputFile)
    i = 0
    j = 0
    fout.write("Column #{0}\n".format(i))
    for line in fin:
        if line[0] != '{' and line[0] != '\n':
            i = line.split()[0]
            j = 0
            fout.write("Scale: {0}\n".format(line.split()[-2].split(":")[1]))
            fout.write("\nColumn #{0}\n".format(i))
        elif line[0] == '{':
            fout.write("Node [{0}]:\t".format(j))
            j += 1
            valuesStr = line.rstrip()[1:-1].split(", ")
            values = map(float, valuesStr)
            maxIndex, maxValue = max(enumerate(values), key=operator.itemgetter(1))
            fout.write("{0} = {1}\n".format(myHash[maxIndex], maxValue))

    fout.close()
    fin.close()

    print "Done!"
    return
#==========================================================================================
'''
Given a file name as a string will attempt to open it for reading.
'''
def openFileWriting(outputFile):
    try:
        fout = open(outputFile,'w')
    except IOError:
        print "Cannot open file to write to: ", outputFile
        return ()

    return fout
#==========================================================================================
'''
Given a file name as a string will attempt to open it for reading.
'''
def openFileReading(outputFile):
    try:
        fout = open(outputFile,'r')
    except IOError:
        print "Cannot open file to read from.."
        return ()

    return fout
#==========================================================================================
'''
Given the sys.argv it will fetch the file handles and open files for reading and writing.
Returns: returns ppa tuple if successful, else the empty tuple.
'''
def parseArgs(args):
    if len(args) != 3:
        print "Usage: <inpuntFile> <outputFile>"
        return ()

    inputFile = args[1]
    outputFile = args[2]

    return inputFile, outputFile
#==========================================================================================
if __name__ == "__main__":
    main()
