F = open("features_list.txt")
f = F.readlines()
fields=[]
# operons = []
# operon_start_sites = []
# operon_end_sites = []
strands = []
locus_tags = []
start_sites = []
end_sites = []
forms = []
features = []


for line in f:
    field = line.strip("\n").split("\t")
    fields.append(field)
    strands.append(field[6])
    start_sites.append(field[3])
    end_sites.append(field[4])
    features.append(field[5])


#-----------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
for i in range(len(fields)):
    for j in range(len(fields)):
        if i != j and start_sites[i] == start_sites[j] and end_sites[i] == end_sites[j] and fields[j][2] != "gene" and fields[j][2] != "pseudogene" and fields[j][2] != "exon":
            fields[i][5] = fields[j][2] #using 2nd index data
            fields[j][5] = "duplicate" #putting data at 5th index and also putiing same data - so fine
exons = [[] for i in range(len(fields))]
for i in range(len(fields)):
    for j in range(len(fields)):
        if i != j and fields[i][2] == "gene" and fields[j][2] == "exon" and start_sites[i] == start_sites[j] and end_sites[i] == end_sites[j]:
            exons[i].append(fields[j][2])
            exons[j].append(fields[i][2])


print("id", "source", "region", "left", "right", "feature", "strand", "new_locus", "old_locus", "exon", sep='\t')
for i in range(len(fields)):
    if len(exons[i]) == 0 and (fields[i][2] == "gene" or fields[i][2] == "pseudogene" or fields[i][2] == "exon" or fields[i][2] == "direct_repeat"):
        print(*fields[i], *exons[i], sep = "\t")
    if len(exons[i]) > 0 and exons[i][0] != "gene" and (fields[i][2] == "gene" or fields[i][2] == "pseudogene"):
        print(*fields[i], *exons[i], sep = "\t")
