#importing relevant packages
from Bio import SeqIO
from Bio import AlignIO
from scipy.interpolate import make_interp_spline, interp1d
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('seaborn-whitegrid')

#################CONSTANT VARIABLES################
PATH = '/Users/kids/Documents/Stanford/Senior Year/Honors Thesis/genetic_alignment_fastas/' #path name of mutliple sequence alignment
ALIGNMENT1a = "dtinc_front_consensus.fasta"
ALIGNMENT1b = "dtinc_revcomp_middle_consensus.fasta"
ALIGNMENT1c = "dtinc_back_consensus.fasta" #name of MSA
ALIGNMENT2a = "osyl_front_consensus.fasta"
ALIGNMENT2b = "osyl_revcomp_back_consensus.fasta"
ALIGNMENT2c = "osyl_revcomp_middle_consensus.fasta" #name of MSA
WINDOW_SIZE = 50 #size of sliding window 




##################FUNCTIONS#######################
#function that returns avg of pairwise hamming distances between sequences
	#input = list of strings
	#will iterate through each sequence and perform pairwise hamming distance calculation
def hamming_distance(stringlist):
	pairwise_hds = [] #list that we will grow w/ pairwise hamming distances
	for i in xrange(len(stringlist)):
		for j in xrange(i+1, len(stringlist)):
			string1 = stringlist[i]
			string2 = stringlist[j]
			pairwise_hds.append(pairwise_hd(string1, string2))
	return sum(pairwise_hds)/len(pairwise_hds) #calculating average of all pairwise hamming distances


#function that returns the Hamming distance between string1 and string2
	#string1 and string2 should be the same length.
def pairwise_hd(string1, string2): 
	distance = 0 #to append as hamming distance is calculated
	for i in range(len(string1)):
		if string1[i] != string2[i]: #if not the same at this site then add to distance counter
			distance += 1
	return distance


#function that puts previous functions together to calculate the variabiity of each alignment 
def calculate_variability(alignment_name):
	#load alignment
	alignment = AlignIO.read(open(PATH + alignment_name), "fasta") #for help, look up fucntion name + documentation
	#make a list of all of our hamming distances to plot later
	hds = []
	#scan through the MSA
	for i in xrange(alignment.get_alignment_length()-WINDOW_SIZE+1): #stop before the very end
		seqslist = []
		for record in alignment:
			#print("%s - %s" % (record.seq[i:i+WINDOW_SIZE], record.id))
			seqslist.append(str(record.seq[i:i+WINDOW_SIZE]))
		#print(seqslist)
		hds.append(hamming_distance(seqslist)) #calculate hamming distance of stored window sequence
	print(hds)


	#makeing a smoothed dataset to plot
	#smoothing function
	#don't use spline - drastically changes the results
	#xcoords = range(len(HDS_all))
	#hds_array = np.array(HDS_all)
	#xcoords_array = np.array(xcoords)

	#x_new = np.linspace(xcoords_array.min(), xcoords_array.max(),150)

	#f = interp1d(xcoords_array, hds_array, kind='quadratic')
	#hds_smooth=f(x_new)

	return (hds)

#######################CODE TO EXECUTE###################
#knit together front, middle, and back of first alignment 
hds_alignment1 = []
hds_alignment1a = calculate_variability(ALIGNMENT1a)
hds_alignment1.extend(hds_alignment1a)
hds_alignment1b = calculate_variability(ALIGNMENT1b)
hds_alignment1.extend(hds_alignment1b)
hds_alignment1c = calculate_variability(ALIGNMENT1c)
hds_alignment1.extend(hds_alignment1c)

#knit together front, middle, and back of second alignment 
hds_alignment2 = []
hds_alignment2a = calculate_variability(ALIGNMENT2a)
hds_alignment2.extend(hds_alignment2a)
hds_alignment2b = calculate_variability(ALIGNMENT2b)
hds_alignment2.extend(hds_alignment2b)
hds_alignment2c = calculate_variability(ALIGNMENT2c)
hds_alignment2.extend(hds_alignment2c)

 #hds_alignment1a + hds_alignment1b + hds_alignment1c
#hds_alignment2 = calculate_variability(ALIGNMENT2a) + calculate_variability(ALIGNMENT2b) + calculate_variability(ALIGNMENT2c)



plt.plot(hds_alignment1, label='D. tinctorius', color = 'blue', linewidth = 2)
plt.plot(hds_alignment2, label='O. sylvatica', color = 'red', linewidth = 2)
plt.title("Sequence variability accross Saxiphilin (window size: " + str(WINDOW_SIZE) + ")")
plt.ylabel("Average Hamming distance of window")
plt.xlabel("Position in sequence")
plt.legend()
plt.show()


    