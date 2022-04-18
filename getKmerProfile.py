import numpy as np

def kmers(dna, K):
    ## Generates kmer counts for given dna sequence and K
    k_dna = ""
    SIZE = np.power(4, K)
    n = len(dna)
    kmer_profile = np.zeros(SIZE)
    for i in range(0, n - K + 1):
        kmer = dna[i:i + K]
        try:
            kmer_profile[kmer2idx[kmer]] += 1
        except:
            pass
    return kmer_profile

def generate_kmer(k):
  ## Generates all possible k-mers
  if k < 1:
    return []
  if k == 1:
    return ['a', 'c', 'g', 't']
  else:
    l = generate_kmer(k-1)
    new_l = []
    for base in ['a', 'c', 'g', 't']:
      for i in l:
        new_l.append(base + i)
    return new_l

def get_probs(dna):
    n = len(dna)
    return {'a': dna.count('a')/n, 'c': dna.count('c')/n, 't': dna.count('t')/n, 'g': dna.count('g')/n}

def merge_probs(p1, p2):
    return {'a': (p1['a'] + p2['a'])/2, 'c': (p1['c'] + p2['c'])/2, 't': (p1['t'] + p2['t'])/2, 'g': (p1['g'] + p2['g'])/2}

def get_prob_for_word(prob, word):
    pr = 1.
    for c in word:
        pr = pr*prob[c]
    return pr

def getKmerProfile(dna):
    """
    Return the kmer profile of a given dna sequence.
    INPUT: dna: a string of dna sequence
    """
    k = 6
    allKmers = generate_kmer(k)
    kmer2idx = dict(zip(allKmers, range(len(allKmers))))
    ## Kmer profile for virus (counts of each kmer)
    kmer_profile = kmers(dna, k)
    ## Calculate probabilities for each base
    p_a_test_virus = {'a': 0, 'c': 0, 't': 0, 'g': 0}
    pr = get_probs(dna)
    p_a_test_virus = merge_probs(p_a_test_virus, pr)
    ## Normalize kmer profile with probabilities of each kmer
    n_hat = len(dna) - k + 1
    cur = []
    for w in allKmers:
        cur.append(kmer_profile[kmer2idx[w]] - n_hat*get_prob_for_word(pr, w))
    ## 1 feature is the length of the viral sequence
    cur.append(np.log(len(dna)))
    test_centralised_counts_vir = np.array(cur)
    return test_centralised_counts_vir
