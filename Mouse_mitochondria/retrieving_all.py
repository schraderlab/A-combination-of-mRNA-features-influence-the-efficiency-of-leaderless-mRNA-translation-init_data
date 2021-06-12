# -----------------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------
import re
F = open("features.txt")
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
leaderless = []


for line in f:
    field = line.strip("\n").split("\t")
    fields.append(field)
    strands.append(field[3])
    leaderless.append(field[4])
    start_sites.append(field[1])
    end_sites.append(field[2])
# for i in range(len(fields)):
#     print(*fields[i], sep ="\t")

# -----------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------entire_sequences, 50 bases, start_codon and start sites-----------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------
from Bio.Data import CodonTable
from Bio.Seq import Seq
standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
mito_table = CodonTable.unambiguous_dna_by_name["Vertebrate Mitochondrial"]

F4= open("mouse_mito_fasta_without_header.txt")
# F4= open("smegmatis_fasta_without_header.txt")
f4 = F4.read()
# print(f4[0:11])

# entire_sequences = [[] for i in range(len(fields))]
# entire_sequences[0].append("entire_seq")
#
# len_entire_sequences = [[] for i in range(len(fields))]
# len_entire_sequences[0].append("len_entire_seq")

# start_codons = [[] for i in range(len(fields))]
# start_codons[0].append("start_codon")

fifty_bases = [[] for i in range(len(fields))]
# fifty_bases[0].append("fifty_bases")

start_positions = [[] for i in range(len(fields))]
# start_positions[0].append("start_position")

#------------------------------sequences_for_all_diff_forms---------------------------------------------------------------------------------------------
for i in range(0, len(fields), 1):
    if strands[i] == "+" and leaderless[i] == "Leadered":
        fifty_bases[i].append(f4[int(start_sites[i])-26:int(start_sites[i])+24].upper())
        start_positions[i].append(26)
    if strands[i] == "-" and leaderless[i] == "Leadered":
        fifty_bases[i].append(Seq(f4[int(end_sites[i])-25:int(end_sites[i])+25]).reverse_complement().upper())
        start_positions[i].append(26)


    if strands[i] == "+" and leaderless[i] != "Leadered":
        fifty_bases[i].append(f4[int(start_sites[i])-1:int(start_sites[i])+49].upper())
        start_positions[i].append(1)
    if strands[i] == "-" and leaderless[i] != "Leadered":
        fifty_bases[i].append(Seq(f4[int(end_sites[i])-50:int(end_sites[i])]).reverse_complement().upper())
        start_positions[i].append(1)

# -----------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------final printing----------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------
for i in range(len(fields)):
    print(*fields[i], *fifty_bases[i], *start_positions[i], sep ="\t")
