
import urllib2
import math
class Homology(object):
    
    def __init__(self, Markov):
        #self.seq = seq
        self.Markov = Markov
        self.dictionary = {}
        self.cntr = 0

        
    def getAlphabet(self, m, level, res, dictionary):
        if level == m:
            return []
        temp = []
        for item in res:
            temp += [i+item for i in ["A","G","C","T"]]
        res = temp[:]
        dictionary[level]={}
        for segm in res:
            dictionary[level][segm] = 0
        self.getAlphabet(m, level+1, res, dictionary)
        return
    
        """supres=[]
        getAlphabet(3,0, [""],dictionary)
        
        print dictionary"""
        
    def makeDictionary(self):
        self.getAlphabet(self.Markov, 0, [""], self.dictionary)


##    def updateCounts(self, seq):
##
##        #self.getAlphabet(self.Markov, 0, [""], self.dictionary)
##
##        for i in range(self.Markov):
##            for j in range(len(seq)-i):
##                self.dictionary[i][seq[j:j+i+1]]+=1
    def updateCounts(self, seq):

        #self.getAlphabet(self.Markov, 0, [""], self.dictionary)

        """self.cntr += len(seq)"""
        for j in range(len(seq) - self.Markov):
            for i in range(self.Markov):
                for i in range(5):
                    try:
                        self.dictionary[i][seq[j:j+i+1]]+=1
                        break
                    except KeyError:
                        pass


        #print self.dictionary
    def getDistribution(self, dictionary):

        for i in range(self.Markov):
            self.cntr = 0
            for alph in self.dictionary[i]:
                    self.cntr += self.dictionary[i][alph]
            for alph in self.dictionary[i]:
                self.dictionary[i][alph] = float(self.dictionary[i][alph])/self.cntr
                if self.dictionary[i][alph] == 0:
                    self.dictionary[i][alph] = 0.0000001
    def getLength (self):
        return self.cntr

    def getPrint(self):
        self.getDistribution(self.dictionary)
        return self.dictionary
"""    
a = {}    
test = Homology("GAATGGAGTT",2)


a= test.getDistribution()
test = Homology("AAAAAAAA",2)
b = test.getDistribution()

print a, b"""


def KL(dis1, dis2):
    res = []
    markov = len(dis1)
    kl = 0
    for m in range(markov):
        for item in dis1[m]:
            kl += (dis1[m][item] * math.log(dis1[m][item]/dis2[m][item],2))
        res.append(kl)
    return res
        
        
Markov = 1        
    

"""s="GAATGGAGTTATATGGACTTCCTGGAGCATCTGCTTCATGAAGAAAAACTGGCACGTCATCAACGTAAACAGGCGATGTATACCCGAATGGCAGCCTTCCCGGCGGTGAAAACGTTCGAAGAGTATGACTTCACATTCGCCACCGGAGCACCGCAGAAGCAACTCCAGTCGTTACGCTCACTCAGCTTCATAGAACGTAATGAAAATATCGTATTACTGGGGCCATCAGGTGTGGGGAAAACCCATCTGGCAATAGCGATGGGCTATGAAGCAGTCCGTGCAGGTATCAAAGTTCGCTTCACAACAGCAGCAGATCTGTTACTTCAGTTATCTACGGCACAACGTCAGGGCCGTTATAAAACGACGCTTCAGCGTGGAGTAATGGCCCCCCGCCTGCTCATCATTGATGAAATAGGCTATCTGCCGTTCAGTCAGGAAGAAGCAAAGCTGTTCTTCCAGGTCATCGCTAAACGTTACGAAAAGAGCGCAATGATCCTGACATCCAATCTGCCGTTCGGGCAGTGGGATCAAACGTTCGCCGGTGATGCAGCACTGACCTCAGCGATGCTGGACCGTATCTTACACCACTCACATGTCGTTCAAATCAAAGGAGAAAGCTATCGACTCAGACAGAAACGAAAGGCCGGGGTTATAGCAGAAGCTAATCCTGAGTAAAACGGTGGATCAATATTGGGCCGTTGGTGGAGATATAAGTGGATCACTTTTCATCCGTCGTTGACA"

#region="NM_009417.2"
region = 'NZ_ACNS01000004.1'
data = urllib2.urlopen('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id='+region+'&rettype=fasta&retmode=text')
data.readline()
pre = ""
Markov = 9
Temp = Homology(Markov)
Temp.makeDictionary()
prels = []
for line in data:
    seq = "".join(i for i in prels) + data.readline().strip()
    prels = list(seq[-Markov:])
    Temp.updateCounts(seq)

seqDis = Temp.getPrint()
#print seqDis

Temp2 = Homology(Markov)

Temp2.makeDictionary()
Temp2.updateCounts(sequence)
geneDis = Temp2.getPrint()
#print geneDis
Divls = KL(seqDis, geneDis)
print "length seq: ", Temp.getLength()
print "my calculation: ", math.log(Temp.getLength(),2)/ Divls[-1]
print Divls"""
header = []

    
filename = "/Users/mitralab/Documents/DNA_Project/hgt_seqs/hgt.txt"
outfilename = "/Users/mitralab/Documents/DNA_Project/hgt_seqs/out.txt"
hgt = open(filename)
out = open(outfilename, 'w')
out.write("sequnce length \t my result \n")
for eachline in hgt:
    if eachline[0] == ">":
        header = eachline.strip().split()
        region = header[-1]
    else:
        sequence = eachline.strip()

        try:
            data = urllib2.urlopen('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id='+region+'&rettype=fasta&retmode=text')
        except URLError:
            data = urllib2.urlopen('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id='+region+'&rettype=fasta&retmode=text')

        data.readline()
        pre = ""
        
        Temp = Homology(Markov)
        Temp.makeDictionary()
        prels = []
        for line in data:
            seq = "".join(i for i in prels) + data.readline().strip()
            prels = list(seq[-Markov:])
            Temp.updateCounts(seq)

        seqDis = Temp.getPrint()
        #print seqDis

        Temp2 = Homology(Markov)

        Temp2.makeDictionary()
        Temp2.updateCounts(sequence)
        geneDis = Temp2.getPrint()
        #print geneDis
        Divls = KL(seqDis, geneDis)
        out.write("%.f \t %.f \n" %(Temp2.getLength(),math.log(Temp.getLength(),2)/ Divls[-1]))
        

out.close()

##print "length seq: ", "length seq: "
##print "my calculation: ", math.log(Temp.getLength(),2)/ Divls[-1]
##print Divls
##        
