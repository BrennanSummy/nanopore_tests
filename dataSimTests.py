import numpy as np
import matplotlib.pyplot as plt

#vvvvvvvvvvvvvvvvvvvv NAMING CONVENTION vvvvvvvvvvvvvvvvvvvv
# Functions are named with: oneTwoThree
# Classes   are named with: OneTwoThree
# Variables are named with: one_two_three
#^^^^^^^^^^^^^^^^^^^^ NAMING CONVENTION ^^^^^^^^^^^^^^^^^^^^

#vvvvvvvvvvvvvvvvvvvv . vvvvvvvvvvvvvvvvvvvv
#^^^^^^^^^^^^^^^^^^^^ . ^^^^^^^^^^^^^^^^^^^^



class KmerHandler():
    '''Handles Kmer map files, creating and storing appropriate dictionaries'''
    def __init__(self, map_filepath,k):
        # Read the map file and create a dictionary for k-mer string keys and corresponding Z scores
        self.kmer_dict = {}
        with open(map_filepath,'r') as f:
            while True:
                line = f.readline()
                if not line:
                    break
                kmer_string, z = self.readLine(line,k)
                self.kmer_dict[kmer_string] = z



    def readLine(self,line_string,k):
        kmer_string = line_string[0:k]

        z = line_string[k+1:-1]

        z = float(z)

        return kmer_string, z

    def printDict(self,n0,n1):
        '''
        Print the items of the kmer to Z score dictionary in the specified range
        Args:
            n0: The first item to print
            n1: The last item to print
        '''
        print(list(self.kmer_dict.items())[n0:n1])
    
    def getDict(self):
        return self.kmer_dict

class FauxDataGenerator():
    '''Given a filepath to a sequence, stores that sequence with the ability to generate faux data given a kmer_dict as input'''
    def __init__(self, sequence_filepath, k, kmer_z_dict, global_mean, global_std_dev):
        with open(sequence_filepath, 'r') as f:
            self.sequence = f.read().replace('\n','')
        print("Sequence file read successfully")
        self.k = k
        self.kmer_expected_current_dict = {kmer_string : self.zToValue(z_score,global_mean,global_std_dev)
                                           for kmer_string, z_score in kmer_z_dict.items()}
        self.mean = global_mean
        self.std_dev = global_std_dev
    
    def genData(self,dwell_time):
        n = len(self.sequence)
        kmer_i = ''
        lmer_i = ''
        expected_current_i = 0
        current_i = 0
        current_trace = []

        # for each letter in sequence
        for i,X in enumerate(self.sequence):
            # if the full 9mer is in the sequence
            if(n-i >= self.k):
                # get the 9mer of relevance by slicing the sequence
                kmer_i      = self.sequence[i:i+self.k]
                expected_current_i   = self.kmer_expected_current_dict[kmer_i]
            # else (meaning that the sequence is cut off, leaving only n letters)
            else:
                # take an average over all 9mers that start with the n accessible letters, and return that current
                lmer_i      = self.sequence[i:]
                expected_current_i   = self.compromiseMer(self.kmer_expected_current_dict,self.k,lmer_i)
            for substep in range(dwell_time):
                current_i = np.random.normal(expected_current_i,0.2*self.std_dev)
                current_trace += [current_i]
        return np.array(current_trace)
                

    def compromiseMer(self, kmer_dict, k, available_sequence):
        '''returns the mean expected current of all possible kmers that start with available_sequence'''
        seq = available_sequence
        for i in range(k-len(available_sequence)):
            seq = seq + 'A'
        print(seq)
        mean_z = kmer_dict[seq]
        return mean_z

    def kMerTree(self, available_sequence, k):
        kmer_array = np.full(4**(k-len(available_sequence)), "", dtype=str)
        # Get all possible kmers by looping through the indices with ACTG using modulus

        return kmer_array

    # Z scores are given in the k-mer tables. This function goes from Z score to actual value
    def zToValue(self, z, current_mean, current_std_dev):
        return (z*current_std_dev + current_mean)

    
k=9
myMap = KmerHandler('9mer_levels_v1.txt',k=k)
myMap.printDict(0,10)

dwell_time = 5
dataGen = FauxDataGenerator('c_elegans_chr3.txt',k=k,kmer_z_dict=myMap.getDict(),global_mean=70,global_std_dev=3)
data = dataGen.genData(dwell_time)
sample_rate = 400 # Base pairs per second
step_rate = sample_rate * dwell_time

plt.plot(data)
plt.title("Simulated ONT Data")
plt.xlabel(f"Timestep ({step_rate} steps per second)")
plt.ylabel("pA")
plt.show()