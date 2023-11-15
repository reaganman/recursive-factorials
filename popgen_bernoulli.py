import find_orf
import translate
from Bio import SeqIO
import sys


# class to phenotype and track samples from fasta file
class Sample:
    def __init__(self, record):
        sampleID_and_date = record.id.split('_')
        self.sampleID = sampleID_and_date[0]
        self.date = sampleID_and_date[1]
        self.sequence = record.seq
        self.translation = self.sequence.translate()
        self.pheno = self.get_pheno()

    def get_pheno(self):
        pos_4 = self.translation.upper()[3]
        if pos_4 == "S":
            return "blue"
        if pos_4 == "R":
            return "orange"
        return "not blue or orange!"


#function to calculate n and k for bernouli_trial and write results
def get_stats(samples, freq, out_file):
    n = len(samples)
    orange_samples = []
    for sample in samples: #find all samples with orange phenotype
        if sample.pheno == "orange":
            orange_samples.append(sample.sampleID)
    k = len(orange_samples)
    p = freq
    bern_stat = bernouli_trial(k, n, freq)

    lines_to_write = [f'Results\n\n', f'p (the frequency of "orange" in the population = {p}\n', f'n (the number of sampled individuals) = {n}\n', f'k (the number of "orange" individuals in the sample set) = {k}\n\n', f'Probability of collecting {n} individuals with {k} being "orange" (given a population frequency of {p}) = {bern_stat}']

    with open(out_file, 'w') as fout:
        for line in lines_to_write:
            fout.write(line)






def factorial(n):
    if n == 1:
        return n
    else:
        return n * factorial(n-1)


#k is number of occurances, n is sample size
def bernouli_trial(k, n, p):
    q = 1 - p
    return (factorial(n)/(factorial(n-k)*factorial(k))) * (p**k) * (q**(n-k))



if __name__=="__main__":
    if len(sys.argv) != 4:
        print("usage: popgen_bernoulli.py <input_fasta> <freq> <out_file>")
        sys.exit()
    infile = sys.argv[1]
    freq = float(sys.argv[2])
    outfile = sys.argv[3]
    samples= []
    for record in SeqIO.parse(infile, 'fasta'): #create sample object for each record in fasta
        sample = Sample(record)
        samples.append(sample)
    get_stats(samples, freq, outfile)

