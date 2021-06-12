# # -----------------------------------------------------------------------------------------------------------------------------------------------------------------
# #------------------------------second way----------------------------------------------------------------------------------------------------------------------------
# # -----------------------------------------------------------------------------------------------------------------------------------------------------------------
import re
# F = open("operon_data_smegmatis.txt")
F = open("operon_data_H37Rv.txt")
f = F.readlines()
fields=[]

strands = []
locus_tags = []
start_sites = []
end_sites = []
forms = []


for line in f:
    field = line.strip("\n").split("\t")
    fields.append(field)



fields_2 = []
for i in range(len(fields)):
    for j in range(len(fields[i])):
        if fields[i][j][0:4] != "back":
            fields_2.append(fields[i][j])



directon_coords_indexes = []
for i in range(len(fields_2)):
    if fields_2[i][0:8] == "Directon":
        directon_coords_indexes.append(i)
directon_coords_indexes.append(len(fields_2))


operons = []
for i in range(len(directon_coords_indexes)-1):
    operons.append(fields_2[directon_coords_indexes[i]:directon_coords_indexes[i+1]])

total_number_of_genes = 0
for i in range(len(operons)):
    number_of_genes_in_each_operon = len(operons[i])-1
    # total_number_of_genes = total_number_of_genes + number_of_genes_in_each_operon
    total_number_of_genes += number_of_genes_in_each_operon
# print(total_number_of_genes)

genes = []
# genes = [[]for i in range(len(total_number_of_genes))]
for i in range(len(operons)):
    for j in range(1, len(operons[i]), 1):
        # print(i+1, *re.split(': |- ', operons[i][0]), operons[i][j][:10], operons[i][j][11:], sep = "\t") #the operons[i][j] needs to be changed based on number of characters
        genes.append([i+1, *re.split(': |- ', operons[i][0]), operons[i][j][:7], operons[i][j][7:]]) #this indexing needs to be changed for all the operons based on length of genes characters
# print(genes)
# print(len(genes))

genes_2 = []
for i in range(len(genes)):
    if genes[i] in genes_2:
        pass
    else:
        genes_2.append(genes[i])
    # for j in range(len(genes_2)):
    #     if genes[i][4] == genes_2[j][4] and genes[i][0] == genes_2[j][0]:
    #         pass
    #     else:
    #         genes_2.append(genes[i])
            # print(genes_2)
# print(genes_2)
# print(len(genes_2))

genes_2 = sorted(genes_2, key=lambda gene: gene[4])
# genes_2 = sorted(genes_2, key=lambda gene: gene[0])
# print(genes_2)
# print(len(genes_2))
genes_2_0_length = len(genes_2[0])
# # print(genes_2_0_length)
#
for i in range(len(genes_2)):
    for j in range(len(genes_2)):
        if j != i and genes_2[i][4] == genes_2[j][4] and genes_2[i][0] != genes_2[j][0]:
            if len(genes_2[i]) < genes_2_0_length:
                genes_2[i].append("BAD")
            if len(genes_2[i]) == genes_2_0_length+1:
                genes_2[i][-1] = "BAD"
        elif len(genes_2[i]) < genes_2_0_length+1:
            genes_2[i].append("GOOD")
# print(genes_2)


bad_operon_numbers = []
for i in range(len(genes_2)):
    if genes_2[i][-1] == "BAD" and genes_2[i][0] not in bad_operon_numbers:
        bad_operon_numbers.append(genes_2[i][0])
    else:
        pass
# print(bad_operon_numbers)

bad_operons = []
for i in range(len(genes_2)):
    for j in range(len(bad_operon_numbers)):
        if genes_2[i][0] == bad_operon_numbers[j]:
            bad_operons.append(genes_2[i])
# for i in range(len(bad_operons)):
#     print(*bad_operons[i], sep = "\t")
# print(bad_operons)
# print(len(bad_operons))


# -----------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------
bad_operon_indexes = []
for i in range(len(bad_operons)):
    if bad_operons[i][6] == "BAD":
        bad_operon_indexes.append(i)
# print(bad_operon_indexes)

bad_operon_indexes_2 = []
for i in range(len(bad_operon_indexes)):
    if i == 0:
        bad_operon_indexes_2.append(bad_operon_indexes[i])
    if i > 0 and bad_operon_indexes[i] > bad_operon_indexes[i-1]+1:
        bad_operon_indexes_2.append(bad_operon_indexes[i])
# print(bad_operon_indexes_2)



#-----------------------------------------------------------------------------------------------------------
index_pairs = []
for i in range(len(bad_operon_indexes_2)):
    if i < len(bad_operon_indexes_2)-1:
        index_pairs.append(bad_operon_indexes[bad_operon_indexes.index(bad_operon_indexes_2[i]):bad_operon_indexes.index(bad_operon_indexes_2[i+1])])
    if i == len(bad_operon_indexes_2)-1:
        index_pairs.append(bad_operon_indexes[bad_operon_indexes.index(bad_operon_indexes_2[i]):])
# print(index_pairs)


operon_pairs_0 = index_pairs
# print(operon_pairs[0][1])
# operon_pairs = [[] for i in in index_pairs]
for i, j in enumerate(index_pairs):
    for k, l in enumerate(j):
        operon_pairs_0[i][k] = bad_operons[l][0]
# print(operon_pairs_0)
# print(len(operon_pairs_0))

operon_pairs_1 = []
for i, j in enumerate(operon_pairs_0):
    if i < len(operon_pairs_0)-1 and operon_pairs_0[i][-1] >= operon_pairs_0[i+1][0]:
        operon_pairs_1.append(operon_pairs_0[i]+operon_pairs_0[i+1])
    elif i > 0 and operon_pairs_0[i][0] <= operon_pairs_0[i-1][-1]:
        pass
    else:
        operon_pairs_1.append(operon_pairs_0[i])
# print(operon_pairs_1)
# print(len(operon_pairs_1))

operon_pairs_2 = []
for i, j in enumerate(operon_pairs_1):
    if i < len(operon_pairs_1)-1 and operon_pairs_1[i][-1] >= operon_pairs_1[i+1][0]:
        operon_pairs_2.append(operon_pairs_1[i]+operon_pairs_1[i+1])
    elif i > 0 and operon_pairs_1[i][0] <= operon_pairs_1[i-1][-1]:
        pass
    else:
        operon_pairs_2.append(operon_pairs_1[i])
# print(operon_pairs_2)
# print(len(operon_pairs_2))


operon_pairs_3 = []
for i, j in enumerate(operon_pairs_2):
    if i < len(operon_pairs_2)-1 and operon_pairs_2[i][-1] >= operon_pairs_2[i+1][0]:
        operon_pairs_3.append(operon_pairs_2[i]+operon_pairs_2[i+1])
    elif i > 0 and operon_pairs_2[i][0] <= operon_pairs_2[i-1][-1]:
        pass
    else:
        operon_pairs_3.append(operon_pairs_2[i])
# print(operon_pairs_3)
# print(len(operon_pairs_3))


operon_pairs = [[] for i in range(len(operon_pairs_3))]
for i, j in enumerate(operon_pairs_3):
    for k, l in enumerate(j):
        if operon_pairs_3[i][k] not in operon_pairs[i]:
            operon_pairs[i].append(operon_pairs_3[i][k])
# print(operon_pairs)
# print(len(operon_pairs))

#-----------------------------------------------------------------------------------------------------------
bad_operons = sorted(bad_operons, key=lambda operon: operon[0])
# for i in range(len(bad_operons)):
#     print(*bad_operons[i], sep = "\t")
bad_operon_first_indexes = []
for i, j in enumerate(operon_pairs):
    for k, l in enumerate(j):
        for m in range(len(bad_operons)):
            if k == 0 and bad_operons[m][0] == operon_pairs[i][k]:
                bad_operon_first_indexes.append(m)
# print(bad_operon_first_indexes)

bad_operon_first_indexes_2 = []
bad_operon_first_indexes_2.append(bad_operon_first_indexes[0])
for i in range(1, len(bad_operon_first_indexes)):
    if bad_operon_first_indexes[i] != bad_operon_first_indexes[i-1]+1:
        bad_operon_first_indexes_2.append(bad_operon_first_indexes[i])
    else:
        pass
# print(bad_operon_first_indexes_2)
# print(len(bad_operon_first_indexes_2))

bad_operon_last_indexes = []
for i in range(1, len(bad_operon_first_indexes_2)):
    bad_operon_last_indexes.append(bad_operon_first_indexes_2[i]-1)
bad_operon_last_indexes.append(len(bad_operons)-1)
# print(bad_operon_last_indexes)
# print(len(bad_operon_last_indexes))
#-----------------------------------------------------------------------------------------------------------
for i, j in enumerate(operon_pairs):
    for k, l in enumerate(j):
        for m in range(len(bad_operons)):
            if k != 0 and bad_operons[m][0] == operon_pairs[i][k]:
                # n = bad_operons.index()
                # bad_operons[m][2] = bad_operons[bad_operons.index(operon_pairs[i][0])][2]
                bad_operons[m][0] = operon_pairs[i][0]
# for i in range(len(bad_operons)):
#     print(*bad_operons[i], sep = "\t")
#
# #-----------------------------------------------------------------------------------------------------------
for i, j in enumerate(operon_pairs):
    # for k, l in enumerate(j):
    for m in range(len(bad_operons)):
        if bad_operons[m][0] == operon_pairs[i][0]:
            bad_operons[m][2] = bad_operons[bad_operon_first_indexes_2[i]][2]
            bad_operons[m][3] = bad_operons[bad_operon_last_indexes[i]][3]

bad_operons_no_dups = []
for i in range(len(bad_operons)):
    if bad_operons[i] not in bad_operons_no_dups:
        bad_operons_no_dups.append(bad_operons[i])

for i in range(len(bad_operons_no_dups)):
    print(*bad_operons_no_dups[i], sep = "\t")

#---------------------------------------------------------------------------------------------------------------
#------------------------------final printing-------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------
for i in range(len(genes_2)):
    for m in range(len(bad_operon_numbers)):
        if genes_2[i][0] == bad_operon_numbers[m]:
            genes_2[i].append("BAD_operon")


for i in range(len(genes_2)):
    if genes_2[i][-1] != "BAD_operon":
        print(*genes_2[i], sep = "\t")
