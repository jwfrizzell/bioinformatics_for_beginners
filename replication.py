# Input:  A string Text and an integer k
# Output: CountDict(Text, k)
def CountDict(Text, k):
	Count = {}
	for i in range(len(Text)-k+1):
		Pattern = Text[i:i+k]
		Count[i] = PatternCount(Pattern, Text)
		#print(str(Pattern) + "\n")
		#print(str(Count[i]) + "\n")

	#print("Exit Count Dict\n\n")
	return Count

# Input:  Strings Pattern and Text
# Output: The number of times Pattern appears in Text
def PatternCount(Pattern, Text):
	count = 0
	for i in range(len(Text)-len(Pattern)+1):
		if Text[i:i+len(Pattern)] == Pattern:
			count = count+1
	return count 

#Get all repeated entries from text based on k-mer length.
#Output new list.
def FrequentWords(Text, k):
	FrequentPatterns = []
	Count = CountDict(Text,k)
	m = max(Count.values())

	for i in Count:
		if Count[i] == m:
			FrequentPatterns.append(Text[i:i+k])
	return remove_duplicates(FrequentPatterns)

#Remove duplicates from list. 
#Output new list with only unique values
def remove_duplicates(list):
	new_list = []
	for item in list:
		if item not in new_list:
			new_list.append(item)
	return new_list

#Reverse the order in which the list has been read. 
#Output string
def reverse(text):
	reversed_list = []
	for i in range(len(text),0, -1):
		reversed_list.append(text[i-1])
	return ''.join(reversed_list)

#Get complements genomes.
#Return the compliments for genomes.
def ReverseComplement(text):
	compliments = []
	for i in range(0,len(text)):
		nucleotide = complement(text[i])
		compliments.append(nucleotide)
	return reverse(''.join(compliments))

# Input:  A character Nucleotide
# Output: The complement of Nucleotide
def complement(nucleotide):
	comp = {
		'A': 'T',
		'a': 't',
		'T': 'A',
		't': 'a',
		'C': 'G',
		'c': 'g',
		'G': 'C',
		'g': 'c'
	}
	return comp.get(nucleotide,'ERROR')

def PatternMatching(Pattern, Genome):
	index = []
	for i in range(0,len(Genome) - len(Pattern)+1):
		if Genome[i:i+len(Pattern)] == Pattern:
			index.append(i)
	return index

def SymbolArray(Genome, symbol):
	array = {}
	n = len(Genome)
	extended_genome = Genome + Genome[0:n//2]
	print(extended_genome)
	for i in range(n):
		array[i] = PatternCount(symbol, extended_genome[i:i+(n//2)])

	return array

def FasterSymbolArray(Genome, symbol):
	array = {}
	n = len(Genome)
	ending_width = n//2
	extended_genome = Genome + Genome[0:ending_width]
	array[0] = PatternCount(symbol, Genome[0:ending_width])
	for i in range(1,n):
		array[i] = array[i-1]
		if extended_genome[i-1] == symbol:
			array[i] = array[i] -1
		if extended_genome[i + ending_width -1] == symbol:
			array[i] = array[i] + 1
	return array

def Skew(Genome):
	skew = {}
	skew[0] = 0
	for i in range(1,len(Genome) +1):
		if Genome[i-1] == 'G':
			skew[i] = skew[i-1] + 1
		elif Genome[i-1] == 'C':
			skew[i] = skew[i-1] -1
		else:
			skew[i] = skew[i-1]
	return skew

def MinimumSkew(Genome):
	positions = []
	skew = Skew(Genome)
	positions = [key for key,value in skew.items() if value == min(skew.values())]
	return positions

def MaximumSkew(Genome):
	positions = []
	skew = Skew(Genome)
	print(skew)
	positions = [key for key,value in skew.items() if value == max(skew.values())]
	return positions

# Get the number of differences between two strings.
def HammingDistance(p, q):
	count = sum(1 for string_one, string_two in zip(p, q) if string_one != string_two)
	return count;

#Get the index position of Patters that are an approximate match
#of the text in sliding window. 
def ApproximatePatternMatching(Pattern, Text, d):
	positions = [] # initializing list of positions
	for i in range(0,len(Text) - len(Pattern)+1):
		if HammingDistance(Pattern, Text[i:i+len(Pattern)]) <= d:
			positions.append(i)
	return positions

# Input:  Strings Pattern and Text, and an integer d
# Output: The number of times Pattern appears in Text with at most d mismatches
def ApproximatePatternCount(Pattern, Text, d):
    count = 0 # initialize count variable
    for i in range(0,len(Text) - len(Pattern)+1):
    	dist = Text[i:i+len(Pattern)]
    	if (HammingDistance(Pattern,Text[i:i+len(Pattern)])) <=d:
    		count = count + 1
    return count

Text="CTTGAAGTGGACCTCTAGTTCCTCTACAAAGAACAGGTTGACCTGTCGCGAAG"
Text2="ATGCCTTACCTAGATGCAATGACGGACGTATTCCTTTTGCCTCAACGGCTCCT"
print(HammingDistance(Text,Text2))
