#-------------------------for this make sure TSS negative from higher to lower-------------
#-------------------------for minus start is actual start and not left for both genes and operons--------------------
#-------------------------make sure "operon_list_smeg.txt" file sorted in order first by genes start and then by operons--------------------
from itertools import repeat

F = open("operon_list_smeg.txt")
f = F.readlines()
fields=[]
operons = []
operon_start_sites = []
operon_end_sites = []
strands = []
locus_tags = []
start_sites = []
end_sites = []
forms = []


for line in f:
    field = line.strip("\n").split("\t")
    fields.append(field)
    operons.append(field[0])
    operon_start_sites.append(field[5])
    operon_end_sites.append(field[6])
    strands.append(field[2])
    locus_tags.append(field[7])
    start_sites.append(field[11])
    end_sites.append(field[12])
    forms.append(field[14])
# print(fields)
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------fist gene_in operon----------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------
if fields[1][2] == "+":
    fields[1][13] = "first_gene"
else:
    fields[1][13] = "in_operon"
# fields[1][13] = "first_gene" #only if positive

for i in range(2, len(fields)):
    if fields[i][2] == "+" and fields[i][0] != "#N/A" and int(fields[i][0]) > int(fields[i-1][0]):
        fields[i][13] = "first_gene"
    elif fields[i][2] == "-" and fields[i][0] != "#N/A" and fields[i+1][0] != "#N/A" and int(fields[i][0]) < int(fields[i+1][0]):
        fields[i][13] = "first_gene"
    elif fields[i][0] == "#N/A" :
        pass
    else:
        fields[i][13] = "in_operon"

# -----------------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------for TRUE----------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------
F2 = open("_6_TSS_plus_smeg")
f2 = F2.readlines()
fields2=[]
for line in f2:
    field2 = line.strip("\n").split()
    fields2.append(field2)
    a = field2[0]
    b = field2[1]


F3 = open("_6_TSS_minus_smeg")
f3 = F3.readlines()
fields3=[]
for line in f3:
    field3 = line.strip("\n").split()
    fields3.append(field3)
    a = field3[0]
    b = field3[1]

# -----------------------------------------------------------------------------------------------------------------------------------------------------------
true = [[] for i in range(len(fields))]
for j in  range (0, len(fields2), 1):
    for i in range (1, len(fields), 1):
        if strands[i] == "+" and fields[i][0] != "#N/A" and float(fields2[j][0]) <= float(operon_start_sites[i]) and float(fields2[j][0]) >= float(operon_start_sites[i]) - 300:
            true[i].append(fields2[j][0]) # this will skip 0 and will directly add in 1, but length of fields and fields3 and fields4 doesnt change
for j in  range (0, len(fields3), 1):
    for i in range (1, len(fields), 1):
        if strands[i] == "-" and fields[i][0] != "#N/A" and float(fields3[j][0]) >= float(operon_end_sites[i]) and float(fields3[j][0]) <= float(operon_end_sites[i]) + 300: # for negative start coordinate is higher and stop is lower so going in direction from left to right-the way wanted in the final version, but noe changed again- so always check.
            true[i].append(fields3[j][0]) # this will skip 0 and will directly add in 1, but length of fields and fields3 and fields4 doesnt change
# print(true)

true_lengths = []
for i in range (0, len(true), 1):
    true_lengths.append(len(true[i]))
# print(lengths3)
# print(max(lengths3))

for i in range (0, len(true), 1):
    if len(true[i]) < max(true_lengths):
        for j in range(max(true_lengths)-len(true[i])):
            true[i].append("x")
# print(true)

# -----------------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------distances---------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------
distances = [[] for i in range(len(fields))]
for i in range(max(true_lengths)):
    distances[0].append("d")

for i in  range (1, len(true), 1):
    for j in range (0, len(true[i]), 1):
        if strands[i] == "+" and true[i][j] != "x":
            distances[i].append(int(start_sites[i])-int(true[i][j]))
        if strands[i] == "+" and true[i][j] == "x":
            distances[i].append("x")

        if strands[i] == "-" and true[i][j] != "x":
            distances[i].append(int(true[i][j])-int(start_sites[i]))
        if strands[i] == "-" and true[i][j] == "x":
            distances[i].append("x")


# -----------------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------diff_forms--------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------
diff_forms = []
for i in  range (1, len(distances), 1):
    for j in range (1, len(distances[i]), 1):
        # if strands[i] == "+" and distances[i][j] != "x" and distances[i][j] <= 25:
        #     diff_forms.append(fields[i][:11] + list(str(j+1)) + fields[i][12:] + true[i] + distances[i]) #list("form_"j) +
        if distances[i][j] != "x" and distances[i][j] <= 25:
            fields.append(fields[i][:14] + list(str(j+1)) + fields[i][15:])
            true.append(true[i])
            distances.append(distances[i])
            strands.append(fields[i][2])
            locus_tags.append(fields[i][7])
            start_sites.append(fields[i][11])
            end_sites.append(fields[i][12])
            forms.append(j+1)




# -----------------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------other_info----------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------
TSS_types = [[] for i in range(len(fields))]
TSS_types[0].append("TSS_type")
numbers_of_TSS = [[] for i in range(len(fields))]
numbers_of_TSS[0].append("number_of_TSS")
internals = [[] for i in range(len(fields))]
internals[0].append("internal")
flags = [[] for i in range(len(fields))]
flags[0].append("flag")
non_coding_RNAs = [[] for i in range(len(fields))]
non_coding_RNAs[0].append("non_coding_RNA")
RNAs = [[] for i in range(len(fields))]
RNAs[0].append("RNA")

for i in range(1, len(fields), 1):
    TSS_types[i].append("NA")
    numbers_of_TSS[i].append("NA")
    internals[i].append("NA")
    flags[i].append("NA")
    non_coding_RNAs[i].append("NA")
    RNAs[i].append("NA")

# -----------------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------entire_sequences, 50 bases, start_codon and start sites-----------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------
from Bio.Data import CodonTable
from Bio.Seq import Seq
standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
mito_table = CodonTable.unambiguous_dna_by_name["Vertebrate Mitochondrial"]

F4= open("smegmatis_fasta_without_header.txt")
f4 = F4.read()


entire_sequences = [[] for i in range(len(fields))]
entire_sequences[0].append("entire_seq")

len_entire_sequences = [[] for i in range(len(fields))]
len_entire_sequences[0].append("len_entire_seq")

start_codons = [[] for i in range(len(fields))]
start_codons[0].append("start_codon")

fifty_bases = [[] for i in range(len(fields))]
fifty_bases[0].append("fifty_bases")

start_positions = [[] for i in range(len(fields))]
start_positions[0].append("start_position")

fifty_bases_genome = [[] for i in range(len(fields))]
fifty_bases_genome[0].append("fifty_bases_genome")

fifty_bases_genome_start_positions = [[] for i in range(len(fields))]
fifty_bases_genome_start_positions[0].append("fifty_bases_genome_start_position")

fifty_bases_genome_2 = [[] for i in range(len(fields))]
fifty_bases_genome_2[0].append("fifty_bases_genome")

fifty_bases_genome_start_positions_2 = [[] for i in range(len(fields))]
fifty_bases_genome_start_positions_2[0].append("fifty_bases_genome_start_position")


#------------------------------based on operon--------------------------------------------------------------------------------------------

for i in range(1, len(fields), 1):
    for j in range (0, max(true_lengths), 1):
        if strands[i] == "+" and true[i][j] != "x" and int(forms[i]) == j+1:
            entire_sequences[i].append(f4[int(true[i][j])-1:int(end_sites[i])].upper())
            len_entire_sequences[i].append(len(entire_sequences[i][0]))
            start_codons[i].append(f4[int(start_sites[i])-1:int(start_sites[i])+2].upper())
            if distances[i][j] >= 25 and len(entire_sequences[i][0]) >= 50:
                fifty_bases[i].append(f4[int(start_sites[i])-26:int(start_sites[i])+24].upper())
                start_positions[i].append(26)
            # if distances[i][j] >= 25 and len(entire_sequences[i][0]) < 50:
            #     fifty_bases[i].append(entire_sequences[i][0])
            #     start_positions[i].append(distances[i][j]+1)
            if distances[i][j] < 25 and len(entire_sequences[i][0]) >= 50:
                fifty_bases[i].append(f4[int(true[i][j])-1:int(true[i][j])+49].upper())
                start_positions[i].append(distances[i][j]+1)
            # if distances[i][j] < 25 and len(entire_sequences[i][0]) < 50:
            #     fifty_bases[i].append(entire_sequences[i][0])
            #     start_positions[i].append(distances[i][j]+1)
            if len(entire_sequences[i][0]) < 50:
                fifty_bases[i].append(entire_sequences[i][0])
                start_positions[i].append(distances[i][j]+1)

        elif strands[i] == "-" and true[i][j] != "x" and int(forms[i]) == j+1:
            entire_sequences[i].append(Seq(f4[int(end_sites[i])-1:int(true[i][j])]).reverse_complement().upper())
            len_entire_sequences[i].append(len(entire_sequences[i][0]))
            start_codons[i].append(Seq(f4[int(start_sites[i])-3:int(start_sites[i])]).reverse_complement().upper())

            if distances[i][j] >= 25 and len(entire_sequences[i][0]) >= 50:
                fifty_bases[i].append(Seq(f4[int(start_sites[i])-25:int(start_sites[i])+25]).reverse_complement().upper())
                start_positions[i].append(26)
            # if distances[i][j] >= 25 and len(entire_sequences[i][0]) < 50:
            #     fifty_bases[i].append(entire_sequences[i][0]).reverse_complement())
            #     start_positions[i].append(distances[i][j]+1)
            if distances[i][j] < 25 and len(entire_sequences[i][0]) >= 50:
                fifty_bases[i].append(Seq(f4[int(true[i][j])-50:int(true[i][j])]).reverse_complement().upper())
                start_positions[i].append(distances[i][j]+1)
            # if distances[i][j] < 25 and len(entire_sequences[i][0]) < 50:
            #     fifty_bases[i].append(entire_sequences[i][0])
            #     start_positions[i].append(distances[i][j]+1)
            if len(entire_sequences[i][0]) < 50:
                fifty_bases[i].append(entire_sequences[i][0])
                start_positions[i].append(distances[i][j]+1)

        else:
            pass

# # -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# #------------------------------for all leaderless RNAs if a TSS not found--------------------------------------------------------------------------------------------------------------------
# # -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
for i in range(1, len(fields), 1):
    if fields[i][15] == "Y" and int(fields[i][14]) == 1 and true[i][0] == "x":
        if strands[i] == "+":
            entire_sequences[i].append(f4[int(start_sites[i])-1:int(end_sites[i])].upper())
            len_entire_sequences[i].append(len(entire_sequences[i][0]))
            start_codons[i].append(f4[int(start_sites[i])-1:int(start_sites[i])+2].upper())
            fifty_bases[i].append(f4[int(start_sites[i])-1:int(start_sites[i])+49].upper())
            start_positions[i].append(1)
            fields[i][16] = "noTSS"
        elif strands[i] == "-":
            entire_sequences[i].append(Seq(f4[int(end_sites[i])-1:int(start_sites[i])]).reverse_complement().upper())
            # print(end_sites[i], start_sites[i])
            len_entire_sequences[i].append(len(entire_sequences[i][0]))
            start_codons[i].append(Seq(f4[int(start_sites[i])-3:int(start_sites[i])]).reverse_complement().upper())
            fifty_bases[i].append(Seq(f4[int(start_sites[i])-50:int(start_sites[i])]).reverse_complement().upper())
            start_positions[i].append(1)
            fields[i][16] = "noTSS"

#changed later after i realized that leaderless(1 nucleotide) not giving proper seq and start position - might not have updated excel file
for i in range(1, len(fields), 1):
    if fields[i][15] == "Y (1 nt leader)" and int(fields[i][14]) == 1 and true[i][0] == "x":
        if strands[i] == "+":
            entire_sequences[i].append(f4[int(start_sites[i])-2:int(end_sites[i])].upper())
            len_entire_sequences[i].append(len(entire_sequences[i][0]))
            start_codons[i].append(f4[int(start_sites[i])-1:int(start_sites[i])+2].upper())
            fifty_bases[i].append(f4[int(start_sites[i])-2:int(start_sites[i])+48].upper())
            start_positions[i].append(2)
            fields[i][16] = "noTSS"
        elif strands[i] == "-":
            entire_sequences[i].append(Seq(f4[int(end_sites[i])-1:int(start_sites[i])+1]).reverse_complement().upper())
            # print(end_sites[i], start_sites[i])
            len_entire_sequences[i].append(len(entire_sequences[i][0]))
            start_codons[i].append(Seq(f4[int(start_sites[i])-3:int(start_sites[i])]).reverse_complement().upper())
            fifty_bases[i].append(Seq(f4[int(start_sites[i])-49:int(start_sites[i])+1]).reverse_complement().upper())
            start_positions[i].append(2)
            fields[i][16] = "noTSS"
# # -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# #------------------------------for all other genes/RNAs--------------------------------------------------------------------------------------------------------------------
# # -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

for i in range(1, len(fields), 1):
    if fields[i][15] != "Y" and fields[i][15] != "Y (1 nt leader)" and int(fields[i][14]) == 1 and true[i][0] == "x":
        if strands[i] == "+":
            entire_sequences[i].append("x")
            len_entire_sequences[i].append("x")
            start_codons[i].append(f4[int(start_sites[i])-1:int(start_sites[i])+2].upper())
            fifty_bases[i].append("x")
            start_positions[i].append("x")
            # fields[i][16] =
        elif strands[i] == "-":
            entire_sequences[i].append("x")
            # print(end_sites[i], start_sites[i])
            len_entire_sequences[i].append("x")
            start_codons[i].append(Seq(f4[int(start_sites[i])-3:int(start_sites[i])]).reverse_complement().upper())
            fifty_bases[i].append("x")
            start_positions[i].append("x")
            # fields[i][16] =

#--------------------------------------------------------------------------------------------------------------------------
for i in range(1, len(fields), 1):
    if len(entire_sequences[i]) <1:
        entire_sequences[i].append("x")
    if len(len_entire_sequences[i]) <1:
        len_entire_sequences[i].append("x")
    if len(start_codons[i]) <1:
        start_codons[i].append("x")
    if len(fifty_bases[i]) <1:
        fifty_bases[i].append("x")
    if len(start_positions[i]) <1:
        start_positions[i].append("x")

# print(len(entire_sequences))
# print(fields)




# -----------------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------combining all----------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------
all = []
all_2 = []
for i in range(len(fields)):
    # print(*fields[i], *true[i], *probable[i], *true_internals[i], *probable_internals[i], sep ="\t") #all the lists are added
    # print(*fields4[i], sep ="\t")
    # all.append([*fields[i], *true[i], *TSS_types[i], *numbers_of_TSS[i], *internals[i], *flags[i], *non_coding_RNAs[i], *RNAs[i],
    #       *distances[i], *entire_sequences[i], *len_entire_sequences[i], *start_codons[i],
    #       *fifty_bases[i], *start_positions[i],])
    all.append([*fields[i], *true[i], *TSS_types[i], *numbers_of_TSS[i], *internals[i], *flags[i], *non_coding_RNAs[i], *RNAs[i],
          *distances[i], *len_entire_sequences[i], *start_codons[i],
          *fifty_bases[i], *start_positions[i]])

# -----------------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------primary_multiform_processed_leaderless----------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------

for i in range(len(all)):
    if (all[i][15] == "Y" or all[i][15] == "Y (1 nt leader)") and int(all[i][14]) == 1 and all[i][17] != "x" and (int(all[i][11]) == int(all[i][17]) or int(all[i][11]) == int(all[i][17])+1 or int(all[i][11]) == int(all[i][17])-1):
        all[i][16] = "primary"
    if (all[i][15] == "Y" or all[i][15] == "Y (1 nt leader)") and int(all[i][14]) == 2 and all[i][18] != "x" and (int(all[i][11]) == int(all[i][18]) or int(all[i][11]) == int(all[i][18])+1 or int(all[i][11]) == int(all[i][18])-1):
        all[i][16] = "multiform"
    if (all[i][15] == "Y" or all[i][15] == "Y (1 nt leader)") and int(all[i][14]) == 3 and all[i][19] != "x" and (int(all[i][11]) == int(all[i][19]) or int(all[i][11]) == int(all[i][19])+1 or int(all[i][11]) == int(all[i][19])-1):
        all[i][16] = "multiform"
    if (all[i][15] == "Y" or all[i][15] == "Y (1 nt leader)") and int(all[i][14]) == 4 and all[i][20] != "x" and (int(all[i][11]) == int(all[i][20]) or int(all[i][11]) == int(all[i][20])+1 or int(all[i][11]) == int(all[i][20])-1):
        all[i][16] = "multiform"
    # if ((all[i][15] == "Y" or all[i][15] == "Y (1 nt leader)") and all[i][17] != "x" and (int(all[i][11]) != int(all[i][17]) or int(all[i][11]) != int(all[i][17])+1 or int(all[i][11]) != int(all[i][17])-1)
    #                                                            and all[i][18] != "x" and (int(all[i][11]) != int(all[i][18]) or int(all[i][11]) != int(all[i][18])+1 or int(all[i][11]) != int(all[i][18])-1)
    #                                                            and all[i][19] != "x" and (int(all[i][11]) != int(all[i][19]) or int(all[i][11]) != int(all[i][19])+1 or int(all[i][11]) != int(all[i][19])-1)
    #                                                            and all[i][17] != "x" and (int(all[i][11]) != int(all[i][20]) or int(all[i][11]) != int(all[i][20])+1 or int(all[i][11]) != int(all[i][20])-1)):
    #     all[i][16] = "processed"
    if ((all[i][15] == "Y" or all[i][15] == "Y (1 nt leader)") and int(all[i][14]) == 1 and all[i][20] != "x" and (int(all[i][11]) != int(all[i][17]) and int(all[i][11]) != int(all[i][17])+1 and int(all[i][11]) != int(all[i][17])-1)
                                                                                                              and (int(all[i][11]) != int(all[i][18]) and int(all[i][11]) != int(all[i][18])+1 and int(all[i][11]) != int(all[i][18])-1)
                                                                                                              and (int(all[i][11]) != int(all[i][19]) and int(all[i][11]) != int(all[i][19])+1 and int(all[i][11]) != int(all[i][19])-1)
                                                                                                              and (int(all[i][11]) != int(all[i][20]) and int(all[i][11]) != int(all[i][20])+1 and int(all[i][11]) != int(all[i][20])-1)):
        # all_2.append([*all[i][:16], "processed", *all[i][17:]])
        all.append([*all[i][:16], "processed", *all[i][17:]])
        # all[i][16] = "processed"
        # print(int(all[i][11]), int(all[i][14]), all[i][15], all[i][16], int(all[i][17]), int(all[i][18]), int(all[i][19]), int(all[i][20]))
    # if ((all[i][15] == "Y" or all[i][15] == "Y (1 nt leader)") and int(all[i][14]) == 1 and all[i][20] == "x" and all[i][19] != "x" and ((int(all[i][11]) != int(all[i][17]) and int(all[i][11]) != int(all[i][18]) and int(all[i][11]) != int(all[i][19])))):
    #                                                                                                                                 # or (int(all[i][11]) != int(all[i][17])+1 and int(all[i][11]) != int(all[i][18])+1 and int(all[i][11]) != int(all[i][19])+1)
    #                                                                                                                                 # or (int(all[i][11]) != int(all[i][17])-1 and int(all[i][11]) != int(all[i][18])-1 and int(all[i][11]) != int(all[i][19])-1))):
    #     # all_2.append([*all[i][:16], "processed", *all[i][17:]])
    #     # all[i][16] = "processed"
    #     print(all[i][11], int(all[i][14]), all[i][15], all[i][16], int(all[i][17]), int(all[i][18]), int(all[i][19]), all[i][20], sep = "\t")
    if ((all[i][15] == "Y" or all[i][15] == "Y (1 nt leader)") and int(all[i][14]) == 1 and all[i][20] == "x" and all[i][19] != "x" and (int(all[i][11]) != int(all[i][17]) and int(all[i][11]) != int(all[i][17])+1 and int(all[i][11]) != int(all[i][17])-1)
                                                                                                                                    and (int(all[i][11]) != int(all[i][18]) and int(all[i][11]) != int(all[i][18])+1 and int(all[i][11]) != int(all[i][18])-1)
                                                                                                                                    and (int(all[i][11]) != int(all[i][19]) and int(all[i][11]) != int(all[i][19])+1 and int(all[i][11]) != int(all[i][19])-1)):
        all.append([*all[i][:16], "processed", *all[i][17:]])
        # all[i][16] = "processed"
        # print(all[i][11], int(all[i][14]), all[i][15], all[i][16], int(all[i][17]), int(all[i][18]), int(all[i][19]), all[i][20])
    if ((all[i][15] == "Y" or all[i][15] == "Y (1 nt leader)") and int(all[i][14]) == 1 and all[i][20] == "x" and all[i][19] == "x" and all[i][18] != "x" and (int(all[i][11]) != int(all[i][17]) and int(all[i][11]) != int(all[i][17])+1 and int(all[i][11]) != int(all[i][17])-1)
                                                                                                                                                          and (int(all[i][11]) != int(all[i][18]) and int(all[i][11]) != int(all[i][18])+1 and int(all[i][11]) != int(all[i][18])-1)):
        all.append([*all[i][:16], "processed", *all[i][17:]])
        # all[i][16] = "processed"
    if ((all[i][15] == "Y" or all[i][15] == "Y (1 nt leader)") and int(all[i][14]) == 1 and all[i][20] == "x" and all[i][19] == "x" and all[i][18] == "x" and all[i][17] != "x" and (int(all[i][11]) != int(all[i][17]) and int(all[i][11]) != int(all[i][17])+1 and int(all[i][11]) != int(all[i][17])-1)):
        all.append([*all[i][:16], "processed", *all[i][17:]])
        # all[i][16] = "processed"

# -----------------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------changing sequences for processed----------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------
for i in range(1, len(all), 1):
    if all[i][16] == "processed" and all[i][2] == "+" and all[i][15] == "Y":
        # entire_sequences[i] = f4[int(all[i][11])-1:int(all[i][12])].upper()
        # len_entire_sequences = len(entire_sequences[i][0])

        # all[i][31] = f4[int(all[i][11])-1:int(all[i][12])].upper()
        all[i][37] = len(f4[int(all[i][11])-1:int(all[i][12])].upper()) #indexes need to be changed based on how many true and distances
        all[i][38] = f4[int(all[i][11])-1:int(all[i][11])+2].upper()
        all[i][39] = f4[int(all[i][11])-1:int(all[i][11])+49].upper()
        all[i][40] = 1

    elif all[i][16] == "processed" and all[i][2] == "+" and all[i][15] == "Y (1 nt leader)":
        # entire_sequences[i] = f4[int(all[i][11])-1:int(all[i][12])].upper()
        # len_entire_sequences = len(entire_sequences[i][0])

        # all[i][31] = f4[int(all[i][11])-2:int(all[i][12])].upper()
        all[i][37] = len(f4[int(all[i][11])-2:int(all[i][12])].upper()) #indexes need to be changed based on how many true and distances
        all[i][38] = f4[int(all[i][11])-1:int(all[i][11])+2].upper()
        all[i][39] = f4[int(all[i][11])-2:int(all[i][11])+48].upper()
        all[i][40] = 2

    elif all[i][16] == "processed" and all[i][2] == "-" and all[i][15] == "Y":
        # entire_sequences = Seq(f4[int(all[i][12])-1:int(all[i][11])]).reverse_complement().upper()
        # print(end_sites[i], start_sites[i])
        # len_entire_sequences = len(entire_sequences[i][0])

        # all[i][31] = Seq(f4[int(all[i][12])-1:int(all[i][11])]).reverse_complement().upper()
        all[i][37] = len(Seq(f4[int(all[i][12])-1:int(all[i][11])]).reverse_complement().upper())
        all[i][38] = Seq(f4[int(all[i][11])-3:int(all[i][11])]).reverse_complement().upper()
        all[i][39] = Seq(f4[int(all[i][11])-50:int(all[i][11])]).reverse_complement().upper()
        all[i][40] = 1

    elif all[i][16] == "processed" and all[i][2] == "-" and all[i][15] == "Y (1 nt leader)":
        # entire_sequences = Seq(f4[int(all[i][12])-1:int(all[i][11])]).reverse_complement().upper()
        # print(end_sites[i], start_sites[i])
        # len_entire_sequences = len(entire_sequences[i][0])

        # all[i][31] = Seq(f4[int(all[i][12])-1:int(all[i][11])+1]).reverse_complement().upper()
        all[i][37] = len(Seq(f4[int(all[i][12])-1:int(all[i][11])+1]).reverse_complement().upper())
        all[i][38] = Seq(f4[int(all[i][11])-3:int(all[i][11])]).reverse_complement().upper()
        all[i][39] = Seq(f4[int(all[i][11])-49:int(all[i][11])+1]).reverse_complement().upper()
        all[i][40] = 2


# -----------------------------------------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------5'UTR----------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------


# -----------------------------------------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------changing Y to N----------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------
for i in range(len(all)):
    if all[i][16] != "primary" and all[i][16] != "multiform" and all[i][16] != "noTSS" and all[i][16] != "processed" and (all[i][15] == "Y" or all[i][15] == "Y (1 nt leader)"):
        all[i][15] = "N"



# -----------------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------TE values----------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------
F5= open("TE_Msmeg.txt")
f5 = F5.readlines()
fields5 = []
for line in f5:
    field5 = line.strip("\n").split("\t")
    fields5.append(field5)

all[0].append("TE_values")
for i in range(len(fields5)):
    for j in range(1, len(all), 1):
        if fields5[i][0] == all[j][7]:
            all[j].append(fields5[i][1])

for i in range(1, len(all)):
    if len(all[i]) < len(all[0]):
        all[i].extend(repeat("x",(len(all[0])-len(all[i]))))




# -----------------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------final printing----------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------
for i in range(len(all)):
    print(*all[i], sep = "\t")
