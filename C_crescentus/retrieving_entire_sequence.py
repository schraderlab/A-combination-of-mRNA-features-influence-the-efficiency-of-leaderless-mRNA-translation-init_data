import collections
from collections import defaultdict
import itertools

F1= open("genes_plus.txt")
f1= F1.readlines()
fields1 = []
operons = []
start_sites = []
end_sites = []
ts_sites = []
locus_tags = []

for line1 in f1:
    field1 = line1.strip("\n").split("\t")
    fields1.append(field1)
    operons.append(field1[0])
    locus_tags.append(field1[1])
    start_sites.append(field1[2])
    end_sites.append(field1[3])
    ts_sites.append(field1[9])
    # a1 = field1[0]
    # b1 = field1[1]
    # c1 = field1[2]
    # d1 = field1[3]
    # e1 = field1[4]
    # f1 = field1[5]
    # g1 = field1[6]
    # h1 = field1[7]
    # i1 = field1[8]
    # print (b1, c1)

# print(start_sites [0], start_sites [1])
# print(end_sites [0], end_sites [1])
# print(ts_sites [0], ts_sites [1])


F2= open("na_1000_fasta_without_header_2")
f2 = F2.read()
# print(f2[(int(ts_sites [1])-1):(int(end_sites[1]))])

# print(ts_sites)
# print(len(ts_sites))
# print(type(ts_sites[1]))
# print(ts_sites[1])
# print(int(ts_sites[1]))
# print(ts_sites[len(ts_sites)-3].isnumeric())
#----------------------------------------for entire_Sequence-----------------------------------------------------
entire_sequences = []

# index = range(1, len(f1), 1)
# for x in index:

for x in range(1, len(f1), 1):
    # print(len(f1))
    if ts_sites[x].isnumeric()  == True  and end_sites[x].isnumeric()  == True:  #this script had used for retrieving start sites, so need to change
        entire_sequences.append(f2[int(ts_sites[x])-1:int(end_sites[x])])
    else:
        entire_sequences.append("N/A")

# print(entire_sequences[0])

for x in range(0, len(entire_sequences), 1):
    # print(len(f1)) #cannot use len(f1) in above range since that is one more than len(entire_sequences)
    print(operons[x+1], locus_tags[x+1], entire_sequences[x], sep="\t")

# --------------------------------------------for different sizes----------------------------------------
# diff_sizes_after_AUG = [15, 30, 45, 60]
# sequences_of_diff_sizes = [[] for i in range(len(diff_sizes_after_AUG))]
#
# for i, j in enumerate(diff_sizes_after_AUG):
#     # print(i, j)
#     for x in range(1, len(f1), 1):
#         sequences_of_diff_sizes[i].append(f2[(int(ts_sites[x])-1):(int(start_sites[x])+j-1)]) #start sites is the number in the list in which numbering starting from 1, but then since we are using the sequence read from file in python, indexing begins at 0. therefore, need to subtract 1
# # print(len(sequences_of_diff_sizes[1]))
#
# for x in range(0, len(entire_sequences), 1):
#     # print(len(f1)) #cannot use len(f1) in above range since that is one more than len(entire_sequences)
#     print(operons[x+1], entire_sequences[x], sequences_of_diff_sizes[0][x], sequences_of_diff_sizes[1][x], sequences_of_diff_sizes[2][x], sequences_of_diff_sizes[3][x])

#--------------------------------------------------------------------------------------
#--------------------------for minus strand---------------------------------------------
#---------------------------------------------------------------------------------------
import collections
from collections import defaultdict
import itertools

F1= open("genes_minus.txt")
f1= F1.readlines()
fields1 = []
operons = []
start_sites = []
end_sites = []
ts_sites = []
locus_tags = []

for line1 in f1:
    field1 = line1.split("\t")
    operons.append(field1[0])
    locus_tags.append(field1[1])
    start_sites.append(field1[2])
    end_sites.append(field1[3])
    ts_sites.append(field1[9])

F2= open("na_1000_fasta_without_header_2_complimentary")
f2 = F2.read()

#--------------------------------------for entire_sequence-----------------------------------------------------
entire_sequences = []

# index = range(1, len(f1), 1)
# for x in index:

for x in range(1, len(f1), 1):
    # print(len(f1))
    if ts_sites[x].isnumeric()  == True  and end_sites[x].isnumeric()  == True:
        entire_sequences.append(f2[int(ts_sites[x])-1:int(end_sites[x])])
    else:
        entire_sequences.append("N/A")

# print(entire_sequences[0])

for x in range(0, len(entire_sequences), 1):
    # print(len(f1)) #cannot use len(f1) in above range since that is one more than len(entire_sequences)
    print(operons[x+1], locus_tags[x+1], entire_sequences[x], sep="\t")

#--------------------------------------------for different sizes----------------------------------------
# diff_sizes_after_AUG = [15, 30, 45, 60]
# sequences_of_diff_sizes = [[] for i in range(len(diff_sizes_after_AUG))]
#
# for i, j in enumerate(diff_sizes_after_AUG):
#     # print(i, j)
#     for x in range(1, len(f1), 1):
#         sequences_of_diff_sizes[i].append(f2[(int(ts_sites[x])-1):(int(start_sites[x])+j-1)])
# # print(len(sequences_of_diff_sizes[1]))
#
# for x in range(0, len(entire_sequences), 1):
#     # print(len(f1)) #cannot use len(f1) in above range since that is one more than len(entire_sequences)
#     print(operons[x+1], entire_sequences[x], sequences_of_diff_sizes[0][x], sequences_of_diff_sizes[1][x], sequences_of_diff_sizes[2][x], sequences_of_diff_sizes[3][x])
