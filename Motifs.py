import random

def Count(Motifs):
	"""
		Count the number of instances of ACGT in the ith
		column. 
	"""
	k = len(Motifs[0])
	keys = ['A','C','G','T']
	count = {key: [0] * k for key in keys}

	k = len(Motifs[0])
	t = len(Motifs)
	for i in range(t):
		for j in range(k):
			symbol = Motifs[i][j]
			count[symbol][j] += 1
	return count

def CountWithPseudocounts(Motifs):
	"""
		Count the number of instances of ACGT in the ith
		column. 
	"""
	k = len(Motifs[0])
	keys = ['A', 'T', 'G','C']
	count = {key: [1] * k for key in keys}
	k = len(Motifs[0])
	t = len(Motifs)
	for i in range(t):
		for j in range(k):
			symbol = Motifs[i][j]
			count[symbol][j] += 1
	return count

def Profile(Motifs):
	"""
		Divide the ith element in dictionary array by the number
		of rows in Motifs.
	"""
	motifs = Count(Motifs)
	first_key = next(iter(motifs))
	profile = {key: [0] * len(motifs[first_key]) for key in motifs.keys()}
	rows = len(Motifs)
	columns = len(Motifs[0])
	for i in motifs.keys():
		for j in range(0,columns):
			profile[i][j] = motifs[i][j] / rows
	return profile

def ProfileWithPseudocounts(Motifs):
	"""
		Divide the ith element in dictionary array by the number
		of rows in Motifs.
	"""
	motifs =  CountWithPseudocounts(Motifs)
	first_key = next(iter(motifs))
	profile = {key: [0] * len(motifs[first_key]) for key in motifs.keys()}
	nucleotide_summation = 0
	columns = len(Motifs[0])
	for i in motifs.keys():
		nucleotide_summation += motifs[i][0]
	for i in motifs.keys():
		for j in range(0,columns):
			profile[i][j] = motifs[i][j] / nucleotide_summation
	return profile

def Consensus(Motifs):
	"""
		Get the most popular nucleotides in each column of the motif 
	"""
	consensus = ""
	count = Count(Motifs);
	first_key = next(iter(count))
	columns = len(count[first_key])
	for i in range(0,columns):
		max_val = 0
		nucleotide = ''
		for key in count.keys():
			if count[key][i] > max_val:
				max_val = count[key][i]
				nucleotide = key
		consensus += nucleotide;
	return consensus;

def Score(Motifs):
	"""
		First construct Consensus(Motifs) and then sum the number of symbols in the 
		j-th column of Motifs that do not match the symbol in position j of the 
		consensus string
	"""
	consensus = Consensus(Motifs)
	rows = len(Motifs)
	columns = len(Motifs[0])
	score = 0;
	for c in range(0, columns):
		for r in range(0, rows):
			if Motifs[r][c] != consensus[c]:
				score +=1
	return score

def Pr(Text,Profile):
	p = 1
	for i in range(0, len(Text)):
		p *= Profile[Text[i]][i]
	return p


def ProfileMostProbablePattern(Text, k, Profile):
	probability = 0
	k_mer = ''
	for i in range(0,len(Text)-k+1):
		pattern = Text[i:i+k]
		p = Pr(pattern,Profile)
		if i == 0:
			probability = p
			k_mer	= pattern
		elif p > probability:
			probability = p
			k_mer = pattern
	return k_mer
    	
def GreedyMotifSearch(Dna, k, t):
	
	#Set BestMotifs equal to the first k-mer from each string in Dna
	BestMotifs = []
	for i in range(0, t):
		BestMotifs.append(Dna[i][0:k])
	#Iterate over the Dna string starting at position 1 in each row of the 
	#matrix. You will then return the highest probablity for each k_mer.
	n = len(Dna[0])
	for i in range(n-k+1):
		Motifs = []
		Motifs.append(Dna[0][i:i+k])
		for j in range(1, t):
			P = Profile(Motifs[0:j])
			Motifs.append(ProfileMostProbablePattern(Dna[j], k, P))
		#Score the Motifs and BestMotifs. If Motifs is a higher score
		#than BestMotifs set the BestMotifs equal to Motifs. 
		if Score(Motifs) < Score(BestMotifs):
			BestMotifs = Motifs
	return BestMotifs

def GreedyMotifSearchWithPseudocounts(Dna, k, t):
	
	#Set BestMotifs equal to the first k-mer from each string in Dna
	BestMotifs = []
	for i in range(0, t):
		BestMotifs.append(Dna[i][0:k])
	#Iterate over the Dna string starting at position 1 in each row of the 
	#matrix. You will then return the highest probablity for each k_mer.
	n = len(Dna[0])
	for i in range(n-k+1):
		Motifs = []
		Motifs.append(Dna[0][i:i+k])
		for j in range(1, t):
			P = ProfileWithPseudocounts(Motifs[0:j])
			Motifs.append(ProfileMostProbablePattern(Dna[j], k, P))
		#Score the Motifs and BestMotifs. If Motifs is a higher score
		#than BestMotifs set the BestMotifs equal to Motifs. 
		if Score(Motifs) < Score(BestMotifs):
			BestMotifs = Motifs
	return BestMotifs

def EntropyPropability(profile):
	"""
		Get the probability for each column and sum the columns and rows.
		M
		sigma -([p(i) * log2(p(i))])
		i = 1 
	"""
	keys = ['A','C','G','T']
	columns = len(profile['A'])
	entropy = [0] * columns
	for i in range(columns):
		column_entropy = 0
		for key in keys:
			pr = profile[key][i]
			if pr != 0.0:
				column_entropy += -(pr * log2(pr))
		entropy[i] = column_entropy
	return sum(entropy)

def Motifs(Profile, Dna):
	"""
		Create list of k-mers for each row in the DNA strand. 
	"""
	k_mers = []
	k = len(Profile['A'])
	for row in Dna:
		k_mers.append(ProfileMostProbablePattern(row,k,Profile))
	return k_mers

def RandomMotifs(Dna,k,t):
	"""
		Generate random motifs from a collection of Dna strands. 
	"""
	k_mers = []
	for row in Dna:
		start_point = random.randint(0,len(row)-k)
		k_mers.append(row[start_point:start_point+k])
	return k_mers

def RandomizedMotifSearch(Dna, k, t):
	"""
		Generate a collection of motifs and set them to BestMotifs. 
		Calculate the profile with pseudocounts for the motifs. 
		If the Score of the motif is less than the best motifs
		set BestMotifs equal to 'M'.
	"""
	M = RandomMotifs(Dna, k, t)
	BestMotifs = M
	while True:
		Profile = ProfileWithPseudocounts(M)
		M = Motifs(Profile, Dna)
		if Score(M) < Score(BestMotifs):
			BestMotifs = M
		else:
			return BestMotifs
		import pdb; pdb.set_trace()  # breakpoint a6506e24 //
		print(BestMotifs)

### DO NOT MODIFY THE CODE BELOW THIS LINE ###
def RepeatedRandomizedMotifSearch(Dna, k, t):
	BestScore = float('inf')
	BestMotifs = []
	for i in range(1):
		Motifs = RandomizedMotifSearch(Dna, k, t)
		CurrScore = Score(Motifs)
		if CurrScore < BestScore:
			BestScore = CurrScore
			BestMotifs = Motifs
	return BestMotifs
0.22 + 0.54 + 0.58 + 0.36  + 0.3
def  Normalize(Probabilities):
	normalized_list = {key: 0 for key in Probabilities.keys()}
	sum_probabilities = sum(Probabilities.values())
	for key, value in Probabilities.items():
		normalized_list[key] =  value / sum_probabilities
	return normalized_list

def WeightedDie(Probabilities):
	k_mer = ''
	p = random.uniform(0,1)
	start_val = 0
	for key, value in Probabilities.items():
		if p >= start_val and p <= (value + start_val):
			k_mer = key
		start_val += value
	return k_mer

def ProfileGeneratedString(Text, profile, k):
	n = len(Text)
	probabilities = {}
	for i in range(0,n-k+1):
		probabilities[Text[i:i+k]] = Pr(Text[i:i+k], profile)
	probabilities = Normalize(probabilities)
	return WeightedDie(probabilities)


def GibbsSampler(Dna, k, t, N):
	BestMotifs = [] 
	motifs = RandomizedMotifSearch(Dna,k, t)
	BestMotifs = motifs
	for j in range(N):
		i = random.randint(0, t-1)
		profile = ProfileWithPseudocounts(motifs[:i] + motifs[(i+1):])
		motifs[i] =  ProfileGeneratedString(Dna[i],profile, len(motifs[0]))
		if Score(motifs) < Score(BestMotifs):
			BestMotifs = motifs
	return BestMotifs



dna = ['ATGAGGTC','GCCCTAGA','AAATAGAT','TTGTGCTA']
motifs = ['GTC','CCC','ATA','GCT']

print(Motifs(Profile(motifs), dna))
#print(RandomizedMotifSearch(dna, k, t))




