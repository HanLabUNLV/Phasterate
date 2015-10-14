'''
Template for creating new script that follow the same general layout.
'''
#==========================================================================================
import sys
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

    currentNumber = 0
    fout = openFileWriting(outputFile)

    for line in fin:
        c = line.split()
        if currentNumber != int(c[0]):
            fout.write("Column {0}| {1}\t{2}\t{3}\t{4}\t{5}\t\n"\
                       .format(c[0], c[2], c[3],c[4],c[5], c[6], c[7]))
            fout.write("Likelihood\tScale\tScale2\n")
            currentNumber = int(c[0])
        fout.write("{0}\t{1}\t{2}\n".format(c[1], c[8].split(":")[1], c[9].split(":")[1]))

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
