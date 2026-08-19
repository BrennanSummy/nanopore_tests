import numpy as np

# Z scores are given in the k-mer tables. This function goes from Z score to actual value
def Z_to_value(Z,mean,std_dev=1):
    return (Z*std_dev + mean)


class KmerHandler():
    def __init__(self, map_filepath,k):
        # Read the map file and create a dictionary for k-mer string keys and corresponding Z scores
        self.kmer_dict = {}
        with open(map_filepath,'r') as f:
            while True:
                line = f.readline()
                if not line:
                    break
                kmer_string, Z = self.readLine(line,k)
                self.kmer_dict[kmer_string] = Z



    def readLine(self,lineString,k):
        kmer_string = lineString[0:k]

        Z = lineString[k+1:-1]
        #if (not Z.isnumeric()):
        #    print(f'Z: {Z} is not a numeric string')
        #    return 0

        Z = float(Z)

        return kmer_string, Z

    def printDict(self,max_entries):
        print(self.kmer_dict[0:max_entries])


myMap = KmerHandler('9mer_levels_v1.txt',9)
myMap.printDict(10)