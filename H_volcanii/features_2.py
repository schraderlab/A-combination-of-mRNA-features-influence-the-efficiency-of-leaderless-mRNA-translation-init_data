#-----------------------------------------------------------------------------------------------------------
#----------------------------------locus file------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
from itertools import repeat


F = open("locus_tag.txt")
f = F.readlines()
locus=[]
for line in f:
    field = line.strip("\n").split("\t")
    locus.append(field)
#-----------------------------------------------------------------------------------------------------------
#-----------------------------------ribo file----------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
F2 = open("ribo.txt")
f2 = F2.readlines()
ribo=[]
for line2 in f2:
    field2 = line2.strip("\n").split("\t")
    ribo.append(field2)
#-----------------------------------------------------------------------------------------------------------
#------------------------------------leaderless_supp file---------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
F3 = open("leaderless_supp.txt")
f3 = F3.readlines()
supp=[]
for line3 in f3:
    field3 = line3.strip("\n").split("\t")
    supp.append(field3)
#-----------------------------------------------------------------------------------------------------------
#--------------------------------------operon file---------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
F4 = open("operon.txt")
f4 = F4.readlines()
operon=[]
for line4 in f4:
    field4 = line4.strip("\n").split("\t")
    operon.append(field4)

#-----------------------------------------------------------------------------------------------------------
#--------------------------------------sequence_file---------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
chromosome = []
F5 = open("chromosome_fasta.txt")
f5 = F5.readlines()
for line5 in f5:
    field5 = line5.strip("\n")
    chromosome.append(field5)
chromosome.pop(0)
# chromosome = ["".join(chromosome[:])]
chromosome = "".join(chromosome[:])
# print(*chromosome)
# print(chromosome)

pHV1 = []
F6 = open("pHV1_fasta.txt")
f6 = F6.readlines()
for i in range(1, len(f6)):
    pHV1.append(f6[i].strip("\n"))
pHV1 = "".join(pHV1[:])
# print(pHV1)

pHV2 = []
F7 = open("pHV2_fasta.txt")
f7 = F7.readlines()
for i in range(1, len(f7)):
    pHV2.append(f7[i].strip("\n"))
pHV2 = "".join(pHV2[:])
# print(pHV2)

pHV3 = []
F8 = open("pHV3_fasta.txt")
f8 = F8.readlines()
for i in range(1, len(f8)):
    pHV3.append(f8[i].strip("\n"))
pHV3 = "".join(pHV3[:])
# print(pHV3)

pHV4 = []
F9 = open("pHV4_fasta.txt")
f9 = F9.readlines()
for i in range(1, len(f9)):
    pHV4.append(f9[i].strip("\n"))
pHV4 = "".join(pHV4[:])
# print(pHV4)
#-----------------------------------------------------------------------------------------------------------
#--------------------------------------locus_and supp-------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
locus[0].append("Replicon")
locus[0].append("HVO (Correct or downstream feature)")
locus[0].append("Length (AA)")
locus[0].append("FC")
locus[0].append("Protein coding start")
locus[0].append("UTR Length")
locus[0].append("type")
for i in range(1, len(supp)):
    for j in range(1, len(locus)):
        if supp[i][5] == locus[j][8] and len(supp[i][8]) != 0 and len(supp[i][10]) == 0 and len(locus[j])==10:
            # locus[j].append([*supp[i][5:10], "leaderless"])
            locus[j].append(supp[i][0])
            locus[j].append(supp[i][5])
            locus[j].append(supp[i][6])
            locus[j].append(supp[i][7])
            locus[j].append(supp[i][8])
            locus[j].append(supp[i][9])
            locus[j].append("leaderless")
        elif supp[i][5] == locus[j][8] and len(supp[i][8]) != 0 and len(supp[i][10]) == 0 and len(locus[j])>10:
            locus.append([*locus[j][0:10], supp[i][0], *supp[i][5:10], "leaderless"])

        elif supp[i][5] == locus[j][8] and len(supp[i][8]) == 0 and len(supp[i][10]) != 0 and len(locus[j])==10:
            # locus[j].append([*supp[i][5:8], *supp[i][10:12], "leadered"])
            locus[j].append(supp[i][0])
            locus[j].append(supp[i][5])
            locus[j].append(supp[i][6])
            locus[j].append(supp[i][7])
            locus[j].append(supp[i][10])
            locus[j].append(supp[i][11])
            locus[j].append("leadered")
        elif supp[i][5] == locus[j][8] and len(supp[i][8]) == 0 and len(supp[i][10]) != 0 and len(locus[j])>10:
            locus.append([*locus[j][0:10], supp[i][0], *supp[i][5:8], *supp[i][10:12], "leadered"])
        # else:
        #     locus[j].append("x")
        #     locus[j].append("x")
        #     locus[j].append("x")
        #     locus[j].append("x")
        #     locus[j].append("x")
        #     locus[j].append("x")
# for i in range(len(locus)):
#     print(*locus[i], sep = "\t")

# locus_length = len(locus[0])
for i in range(1, len(locus)):
    if len(locus[i]) < len(locus[0]):
        locus[i].extend(repeat("x",(len(locus[0])-len(locus[i]))))
        # locus[i].append("x")*(len(locus[0])-len(locus[i]))

locus_tags =  []
for j in range(1, len(locus)):
    locus_tags.append(locus[j][8])

for i in range (1, len(supp)):
    if supp[i][5] not in locus_tags and len(supp[i][8]) != 0 and len(supp[i][10]) == 0:
        locus.append(["SUPP_FILE", "x", "x", "x", "x", "x", supp[i][2], "x", supp[i][5], "x", supp[i][0], *supp[i][5:10], "leaderless"])
    if supp[i][5] not in locus_tags and len(supp[i][8]) == 0 and len(supp[i][10]) != 0:
        locus.append(["SUPP_FILE", "x", "x", "x", "x", "x", supp[i][2], "x", supp[i][5], "x", supp[i][0], *supp[i][5:8], *supp[i][10:12], "leadered"])



# for i in range(len(locus)):
#     print(*locus[i], sep = "\t")

#-----------------------------------------------------------------------------------------------------------
#--------------------------------------locus-supp and operon------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
locus[0].append("#operonID")
for i in range(1, len(operon)):
    for j in range(1, len(locus)):
        if operon[i][1] == locus[j][8] and locus[j][0] != "SUPP_FILE":
            locus[j].append(operon[i][0])
        elif operon[i][1] == locus[j][8] and locus[j][0] == "SUPP_FILE":
            locus[j].append(operon[i][0])
            locus[j][3] = operon[i][2]
            locus[j][4] = operon[i][3]
        # else:
        #     locus[j].append("x")

for i in range(1, len(locus)):
    if len(locus[i]) < len(locus[0]):
        locus[i].extend(repeat("x",(len(locus[0])-len(locus[i]))))


operon_ids =  []
for j in range(1, len(locus)):
    operon_ids.append(locus[j][17])

for i in range (1, len(operon)):
    if operon[i][0] not in operon_ids:
        locus.append(["OPERON_FILE", "x", "x", operon[i][2], operon[i][3], "x", operon[i][4], "x", operon[i][1], "x", "x", "x", "x", "x", "x", "x", "x", operon[i][0]])
    # if operon[i][0] not in operon_ids:
    #     locus.append(["x", "x", "x", operon[i][2], operon[i][3], "x", operon[i][4], "x", operon[i][1], "x", "x", "x", "x", "x", "x", "x", "x", operon[i][0]])



# for i in range(len(locus)):
#     print(*locus[i], sep = "\t")

#-----------------------------------------------------------------------------------------------------------
#---------------------------first gene/operons--------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
locus = sorted(locus, key=lambda start:start[3])
locus = sorted(locus, key=lambda operon:operon[17])

locus[0].append("first_gene or in_operon")

if locus[1][6] == "+" or locus[1][6] == "Fwd":
    locus[1].append("first_gene")
else:
    locus[1].append("in_operon")

for i in range(2, len(locus), 1):
    if(locus[i][6] == "+" or locus[i][6] == "Fwd") and locus[i][17] != "x" and locus[i][17] != locus[i-1][17]:
        locus[i].append("first_gene")
    elif(locus[i][6] == "-" or locus[i][6] == "Rev") and locus[i][17] != "x" and locus[i][17] != locus[i+1][17]:
        locus[i].append("first_gene")
    elif locus[i][17] == "x":
        pass
    else:
        locus[i].append("in_operon")

for i in range(1, len(locus)):
    if len(locus[i]) < len(locus[0]):
        locus[i].extend(repeat("x",(len(locus[0])-len(locus[i]))))

# for i in range(len(locus)):
#     print(*locus[i], sep = "\t")

#-----------------------------------------------------------------------------------------------------------
#---------------------------leaderless and also operon-processed--------------------------------------------
#-----------------------------------------------------------------------------------------------------------
for i in range(1, len(locus), 1):
    if locus[i][16] == "leaderless" and locus[i][18] == "in_operon":
        locus.append(["PROCESSED", *locus[i][1:15], "x", "x", *locus[i][17:]])




#-----------------------------------------------------------------------------------------------------------
#---------------------------sequences-----------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
from Bio.Data import CodonTable
from Bio.Seq import Seq
standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
mito_table = CodonTable.unambiguous_dna_by_name["Vertebrate Mitochondrial"]

# print(locus[1138][17][:4])
# print(locus[1138][10])
# print(locus[1139])
locus[0].append("fifty_bases")
locus[0].append("start_pos")
locus[0].append("start_codon")
for i in range(1, len(locus)):
    if locus[i][6] == "+" or locus[i][6] == "Fwd":
        if (locus[i][10] == "CHR" or locus[i][17][:3] == "CHR") and locus[i][14] != "x" and locus[i][16] != "x" and int(locus[i][15])<25: #locus[i][14] != "x" is extra info, initially i had written locus[i][3] to see if start is present - but later figured starts are given with supp fle - so can use from there directly so that the ones for whom I could not assign start, I could retrieve the sequences
            locus[i].append(chromosome[int(locus[i][14])-int(locus[i][15])-1:int(locus[i][14])-int(locus[i][15])+49])
            # print(chromosome[int(locus[i][14])-int(locus[i][15])-1:int(locus[i][14])-int(locus[i][15])+49])
            locus[i].append(int(locus[i][15])+1)
            locus[i].append(chromosome[int(locus[i][14])-1:int(locus[i][14])+2])
        if (locus[i][10] == "CHR" or locus[i][17][:3] == "CHR") and locus[i][14] != "x" and locus[i][16] != "x" and int(locus[i][15])>=25:
            locus[i].append(chromosome[int(locus[i][14])-26:int(locus[i][14])+24])
            locus[i].append(26)
            locus[i].append(chromosome[int(locus[i][14])-1:int(locus[i][14])+2])

        if (locus[i][10] == "pHV1" or locus[i][17][:4] == "pHV1") and locus[i][14] != "x" and locus[i][16] != "x" and int(locus[i][15])<25:
            locus[i].append(pHV1[int(locus[i][14])-int(locus[i][15])-1:int(locus[i][14])-int(locus[i][15])+49])
            # print(pHV1[int(locus[i][14])-int(locus[i][15])-1:int(locus[i][14])-int(locus[i][15])+49])
            locus[i].append(int(locus[i][15])+1)
            locus[i].append(pHV1[int(locus[i][14])-1:int(locus[i][14])+2])
        if (locus[i][10] == "pHV1" or locus[i][17][:4] == "pHV1") and locus[i][14] != "x" and locus[i][16] != "x" and int(locus[i][15])>=25:
            locus[i].append(pHV1[int(locus[i][14])-26:int(locus[i][14])+24])
            locus[i].append(26)
            locus[i].append(pHV1[int(locus[i][14])-1:int(locus[i][14])+2])

        if (locus[i][10] == "pHV2" or locus[i][17][:4] == "pHV2") and locus[i][14] != "x" and locus[i][16] != "x" and int(locus[i][15])<25:
            locus[i].append(pHV2[int(locus[i][14])-int(locus[i][15])-1:int(locus[i][14])-int(locus[i][15])+49])
            locus[i].append(int(locus[i][15])+1)
            locus[i].append(pHV2[int(locus[i][14])-1:int(locus[i][14])+2])
        if (locus[i][10] == "pHV2" or locus[i][17][:4] == "pHV2") and locus[i][14] != "x" and locus[i][16] != "x" and int(locus[i][15])>=25:
            locus[i].append(pHV2[int(locus[i][14])-26:int(locus[i][14])+24])
            locus[i].append(26)
            locus[i].append(pHV2[int(locus[i][14])-1:int(locus[i][14])+2])

        if (locus[i][10] == "pHV3" or locus[i][17][:4] == "pHV3") and locus[i][14] != "x" and locus[i][16] != "x" and int(locus[i][15])<25:
            locus[i].append(pHV3[int(locus[i][14])-int(locus[i][15])-1:int(locus[i][14])-int(locus[i][15])+49])
            locus[i].append(int(locus[i][15])+1)
            locus[i].append(pHV3[int(locus[i][14])-1:int(locus[i][14])+2])
        if (locus[i][10] == "pHV3" or locus[i][17][:4] == "pHV3") and locus[i][14] != "x" and locus[i][16] != "x" and int(locus[i][15])>=25:
            locus[i].append(pHV3[int(locus[i][14])-26:int(locus[i][14])+24])
            locus[i].append(26)
            locus[i].append(pHV3[int(locus[i][14])-1:int(locus[i][14])+2])

        if (locus[i][10] == "pHV4" or locus[i][17][:4] == "pHV4") and locus[i][14] != "x" and locus[i][16] != "x" and int(locus[i][15])<25:
            locus[i].append(pHV4[int(locus[i][14])-int(locus[i][15])-1:int(locus[i][14])-int(locus[i][15])+49])
            locus[i].append(int(locus[i][15])+1)
            locus[i].append(pHV4[int(locus[i][14])-1:int(locus[i][14])+2])
        if (locus[i][10] == "pHV4" or locus[i][17][:4] == "pHV4") and locus[i][14] != "x" and locus[i][16] != "x" and int(locus[i][15])>=25:
            locus[i].append(pHV4[int(locus[i][14])-26:int(locus[i][14])+24])
            locus[i].append(26)
            locus[i].append(pHV4[int(locus[i][14])-1:int(locus[i][14])+2])



#---------------------------for positive operons--------------------------------------------------------------------------------
        if (locus[i][10] == "CHR" or locus[i][17][:3] == "CHR") and locus[i][3] != "x" and locus[i][16] == "x" and locus[i][18] == "in_operon":
            locus[i].append(chromosome[int(locus[i][3])-26:int(locus[i][3])+24])
            locus[i].append(26)
            locus[i].append(chromosome[int(locus[i][3])-1:int(locus[i][3])+2])

        if (locus[i][10] == "pHV1" or locus[i][17][:4] == "pHV1") and locus[i][3] != "x" and locus[i][16] == "x" and locus[i][18] == "in_operon":
            locus[i].append(pHV1[int(locus[i][3])-26:int(locus[i][3])+24])
            locus[i].append(26)
            locus[i].append(pHV1[int(locus[i][3])-1:int(locus[i][3])+2])


        if (locus[i][10] == "pHV2" or locus[i][17][:4] == "pHV2") and locus[i][3] != "x" and locus[i][16] == "x" and locus[i][18] == "in_operon":
            locus[i].append(pHV2[int(locus[i][3])-26:int(locus[i][3])+24])
            locus[i].append(26)
            locus[i].append(pHV2[int(locus[i][3])-1:int(locus[i][3])+2])


        if (locus[i][10] == "pHV3" or locus[i][17][:4] == "pHV3") and locus[i][3] != "x" and locus[i][16] == "x" and locus[i][18] == "in_operon":
            locus[i].append(pHV3[int(locus[i][3])-26:int(locus[i][3])+24])
            locus[i].append(26)
            locus[i].append(pHV3[int(locus[i][3])-1:int(locus[i][3])+2])


        if (locus[i][10] == "pHV4" or locus[i][17][:4] == "pHV4") and locus[i][3] != "x" and locus[i][16] == "x" and locus[i][18] == "in_operon":
            locus[i].append(pHV4[int(locus[i][3])-26:int(locus[i][3])+24])
            locus[i].append(26)
            locus[i].append(pHV4[int(locus[i][3])-1:int(locus[i][3])+2])

#---------------------------for positive first_genes--------------------------------------------------------------------------------
        if (locus[i][10] == "CHR" or locus[i][17][:3] == "CHR") and locus[i][3] != "x" and locus[i][16] == "x" and locus[i][18] == "first_gene":
            locus[i].append("x")
            locus[i].append("x")
            locus[i].append(chromosome[int(locus[i][3])-1:int(locus[i][3])+2])

        if (locus[i][10] == "pHV1" or locus[i][17][:4] == "pHV1") and locus[i][3] != "x" and locus[i][16] == "x" and locus[i][18] == "first_gene":
            locus[i].append("x")
            locus[i].append("x")
            locus[i].append(pHV1[int(locus[i][3])-1:int(locus[i][3])+2])


        if (locus[i][10] == "pHV2" or locus[i][17][:4] == "pHV2") and locus[i][3] != "x" and locus[i][16] == "x" and locus[i][18] == "first_gene":
            locus[i].append("x")
            locus[i].append("x")
            locus[i].append(pHV2[int(locus[i][3])-1:int(locus[i][3])+2])


        if (locus[i][10] == "pHV3" or locus[i][17][:4] == "pHV3") and locus[i][3] != "x" and locus[i][16] == "x" and locus[i][18] == "first_gene":
            locus[i].append("x")
            locus[i].append("x")
            locus[i].append(pHV3[int(locus[i][3])-1:int(locus[i][3])+2])


        if (locus[i][10] == "pHV4" or locus[i][17][:4] == "pHV4") and locus[i][3] != "x" and locus[i][16] == "x" and locus[i][18] == "first_gene":
            locus[i].append("x")
            locus[i].append("x")
            locus[i].append(pHV4[int(locus[i][3])-1:int(locus[i][3])+2])

#---------------------------for negative--------------------------------------------------------------------------------

    if locus[i][6] == "-" or locus[i][6] == "Rev":
        if (locus[i][10] == "CHR" or locus[i][17][:3] == "CHR") and locus[i][14] != "x" and locus[i][16] != "x" and int(locus[i][15])<25: #as i mentioned before locus[i][14] != "x"-extra
            locus[i].append(Seq(chromosome[int(locus[i][14])+int(locus[i][15])-50:int(locus[i][14])+int(locus[i][15])]).reverse_complement())
            # print(int(locus[i][14])+int(locus[i][15]))
            locus[i].append(int(locus[i][15])+1)
            locus[i].append(Seq(chromosome[int(locus[i][14])-3:int(locus[i][14])]).reverse_complement())
        if (locus[i][10] == "CHR" or locus[i][17][:3] == "CHR") and locus[i][14] != "x" and locus[i][16] != "x" and int(locus[i][15])>=25:
            locus[i].append(Seq(chromosome[int(locus[i][14])-25:int(locus[i][14])+25]).reverse_complement())
            locus[i].append(26)
            locus[i].append(Seq(chromosome[int(locus[i][14])-3:int(locus[i][14])]).reverse_complement())

        if (locus[i][10] == "pHV1" or locus[i][17][:4] == "pHV1") and locus[i][14] != "x" and locus[i][16] != "x" and int(locus[i][15])<25:
            locus[i].append(Seq(pHV1[int(locus[i][14])+int(locus[i][15])-50:int(locus[i][14])+int(locus[i][15])]).reverse_complement())
            locus[i].append(int(locus[i][15])+1)
            locus[i].append(Seq(pHV1[int(locus[i][14])-3:int(locus[i][14])]).reverse_complement())
        if (locus[i][10] == "pHV1" or locus[i][17][:4] == "pHV1") and locus[i][14] != "x" and locus[i][16] != "x" and int(locus[i][15])>=25:
            locus[i].append(Seq(pHV1[int(locus[i][14])-25:int(locus[i][14])+25]).reverse_complement())
            locus[i].append(26)
            locus[i].append(Seq(pHV1[int(locus[i][14])-3:int(locus[i][14])]).reverse_complement())

        if (locus[i][10] == "pHV2" or locus[i][17][:4] == "pHV2") and locus[i][14] != "x" and locus[i][16] != "x" and int(locus[i][15])<25:
            locus[i].append(Seq(pHV2[int(locus[i][14])+int(locus[i][15])-50:int(locus[i][14])+int(locus[i][15])]).reverse_complement())
            locus[i].append(int(locus[i][15])+1)
            locus[i].append(Seq(pHV2[int(locus[i][14])-3:int(locus[i][14])]).reverse_complement())
        if (locus[i][10] == "pHV2" or locus[i][17][:4] == "pHV2") and locus[i][14] != "x" and locus[i][16] != "x" and int(locus[i][15])>=25:
            locus[i].append(Seq(pHV2[int(locus[i][14])-25:int(locus[i][14])+25]).reverse_complement())
            locus[i].append(26)
            locus[i].append(Seq(pHV2[int(locus[i][14])-3:int(locus[i][14])]).reverse_complement())

        if (locus[i][10] == "pHV3" or locus[i][17][:4] == "pHV3") and locus[i][14] != "x" and locus[i][16] != "x" and int(locus[i][15])<25:
            locus[i].append(Seq(pHV3[int(locus[i][14])+int(locus[i][15])-50:int(locus[i][14])+int(locus[i][15])]).reverse_complement())
            locus[i].append(int(locus[i][15])+1)
            locus[i].append(Seq(pHV3[int(locus[i][14])-3:int(locus[i][14])]).reverse_complement())
        if (locus[i][10] == "pHV3" or locus[i][17][:4] == "pHV3") and locus[i][14] != "x" and locus[i][16] != "x" and int(locus[i][15])>=25:
            locus[i].append(Seq(pHV3[int(locus[i][14])-25:int(locus[i][14])+25]).reverse_complement())
            locus[i].append(26)
            locus[i].append(Seq(pHV3[int(locus[i][14])-3:int(locus[i][14])]).reverse_complement())

        if (locus[i][10] == "pHV4" or locus[i][17][:4] == "pHV4") and locus[i][14] != "x" and locus[i][16] != "x" and int(locus[i][15])<25:
            locus[i].append(Seq(pHV4[int(locus[i][14])+int(locus[i][15])-50:int(locus[i][14])+int(locus[i][15])]).reverse_complement())
            locus[i].append(int(locus[i][15])+1)
            locus[i].append(Seq(pHV4[int(locus[i][14])-3:int(locus[i][14])]).reverse_complement())
        if (locus[i][10] == "pHV4" or locus[i][17][:4] == "pHV4") and locus[i][14] != "x" and locus[i][16] != "x" and int(locus[i][15])>=25:
            locus[i].append(Seq(pHV4[int(locus[i][14])-25:int(locus[i][14])+25]).reverse_complement())
            locus[i].append(26)
            locus[i].append(Seq(pHV4[int(locus[i][14])-3:int(locus[i][14])]).reverse_complement())

#---------------------------for negative operons--------------------------------------------------------------------------------

        if (locus[i][10] == "CHR" or locus[i][17][:3] == "CHR") and locus[i][3] != "x" and locus[i][16] == "x" and locus[i][18] == "in_operon":
            locus[i].append(Seq(chromosome[int(locus[i][4])-25:int(locus[i][4])+25]).reverse_complement())
            locus[i].append(26)
            locus[i].append(Seq(chromosome[int(locus[i][4])-3:int(locus[i][4])]).reverse_complement())

        if (locus[i][10] == "pHV1" or locus[i][17][:4] == "pHV1") and locus[i][3] != "x" and locus[i][16] == "x" and locus[i][18] == "in_operon":
            locus[i].append(Seq(pHV1[int(locus[i][4])-25:int(locus[i][4])+25]).reverse_complement())
            locus[i].append(26)
            locus[i].append(Seq(pHV1[int(locus[i][4])-3:int(locus[i][4])]).reverse_complement())

        if (locus[i][10] == "pHV2" or locus[i][17][:4] == "pHV2") and locus[i][3] != "x" and locus[i][16] == "x" and locus[i][18] == "in_operon":
            locus[i].append(Seq(pHV2[int(locus[i][4])-25:int(locus[i][4])+25]).reverse_complement())
            locus[i].append(26)
            locus[i].append(Seq(pHV2[int(locus[i][4])-3:int(locus[i][4])]).reverse_complement())

        if (locus[i][10] == "pHV3" or locus[i][17][:4] == "pHV3") and locus[i][3] != "x" and locus[i][16] == "x" and locus[i][18] == "in_operon":
            locus[i].append(Seq(pHV3[int(locus[i][4])-25:int(locus[i][4])+25]).reverse_complement())
            locus[i].append(26)
            locus[i].append(Seq(pHV3[int(locus[i][4])-3:int(locus[i][4])]).reverse_complement())

        if (locus[i][10] == "pHV4" or locus[i][17][:4] == "pHV4") and locus[i][3] != "x" and locus[i][16] == "x" and locus[i][18] == "in_operon":
            locus[i].append(Seq(pHV4[int(locus[i][4])-25:int(locus[i][4])+25]).reverse_complement())
            locus[i].append(26)
            locus[i].append(Seq(pHV4[int(locus[i][4])-3:int(locus[i][4])]).reverse_complement())


#---------------------------for negative first_genes--------------------------------------------------------------------------------

        if (locus[i][10] == "CHR" or locus[i][17][:3] == "CHR") and locus[i][3] != "x" and locus[i][16] == "x" and locus[i][18] == "first_gene":
            locus[i].append("x")
            locus[i].append("x")
            locus[i].append(Seq(chromosome[int(locus[i][4])-3:int(locus[i][4])]).reverse_complement())

        if (locus[i][10] == "pHV1" or locus[i][17][:4] == "pHV1") and locus[i][3] != "x" and locus[i][16] == "x" and locus[i][18] == "first_gene":
            locus[i].append("x")
            locus[i].append("x")
            locus[i].append(Seq(pHV1[int(locus[i][4])-3:int(locus[i][4])]).reverse_complement())

        if (locus[i][10] == "pHV2" or locus[i][17][:4] == "pHV2") and locus[i][3] != "x" and locus[i][16] == "x" and locus[i][18] == "first_gene":
            locus[i].append("x")
            locus[i].append("x")
            locus[i].append(Seq(pHV2[int(locus[i][4])-3:int(locus[i][4])]).reverse_complement())

        if (locus[i][10] == "pHV3" or locus[i][17][:4] == "pHV3") and locus[i][3] != "x" and locus[i][16] == "x" and locus[i][18] == "first_gene":
            locus[i].append("x")
            locus[i].append("x")
            locus[i].append(Seq(pHV3[int(locus[i][4])-3:int(locus[i][4])]).reverse_complement())

        if (locus[i][10] == "pHV4" or locus[i][17][:4] == "pHV4") and locus[i][3] != "x" and locus[i][16] == "x" and locus[i][18] == "first_gene":
            locus[i].append("x")
            locus[i].append("x")
            locus[i].append(Seq(pHV4[int(locus[i][4])-3:int(locus[i][4])]).reverse_complement())

#---------------------------start codon for all others - not work since no chromosome annotation --------------------------------------------------------------------------------
for i in range(1, len(locus)):
    if len(locus[i]) < len(locus[0]):
        locus[i].extend(repeat("x",(len(locus[0])-len(locus[i]))))

# if (locus[i][6] == "+" or locus[i][6] == "Fwd") and locus[i][21] == "x" and locus[i][21] != "x":
#         # locus[i].append(x)
#         # locus[i].append(26)
#         locus[i] = pHV4[int(locus[i][3])-1:int(locus[i][3])+2]
#
# if (locus[i][6] == "-" or locus[i][6] == "Rev") and locus[i][21] == "x" and locus[i][21] != "x":
#         # locus[i].append(x)
#         # locus[i].append(26)
#         locus[i][21] = Seq(chromosome[int(locus[i][4])-3:int(locus[i][4])]).reverse_complement())
#-----------------------------------------------------------------------------------------------------------
#---------------------------terms---------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
# locus = sorted(locus, key=lambda seq:seq[19])
# print (locus[2][15])
# print(type(locus[2][15]))
# print (locus[3][15])
# print(type(locus[3][15]))
# print (locus[4][15])
# print(type(locus[4][15]))
# print (locus[5][15])
# print(type(locus[5][15]))

locus[0].append("terms")
for i in range(1, len(locus)):
    if locus[i][15] == "0":
        locus[i].append("Leaderless")

    if locus[i][15] != "x" and locus[i][15] != "0" and int(locus[i][15]) != "x" and locus[i][18] == "first_gene":
        locus[i].append("Leadered-first_gene")

    if locus[i][15] != "x" and locus[i][15] != "0" and int(locus[i][15]) != "x" and locus[i][18] == "x":
        locus[i].append("Leadered-not_sure")

    if locus[i][15] != "0" and locus[i][18] == "in_operon":
        locus[i].append("operon")

for i in range(1, len(locus)):
    if len(locus[i]) < len(locus[0]):
        locus[i].extend(repeat("x",(len(locus[0])-len(locus[i]))))

#-----------------------------------------------------------------------------------------------------------
#---------------------------final printing------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------

for i in range(len(locus)):
    print(*locus[i], sep = "\t")




#-----------------------------------------------------------------------------------------------------------
#---------------------------end-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
