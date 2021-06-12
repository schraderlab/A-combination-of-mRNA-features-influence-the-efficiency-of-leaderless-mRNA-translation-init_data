# -----------------------------for plus strand-------------------------
F = open("operon_list_plus.txt")
f = F.readlines()
fields=[]
for line in f:
    field = line.strip("\n").split("\t")
    fields.append(field)
    a = field[0]
    b = field[1]
    c = field[2]
    d = field[3]
    e = field[4]
    f = field[5]
    g = field[6]
    # h = field[7]


F2 = open("_6_TSS_plus")
f2 = F2.readlines()
fields2=[]
for line in f2:
    field2 = line.strip("\n").split()
    fields2.append(field2)
    a = field2[0]
    b = field2[1]
    # c = field2[2]
    # d = field2[3]
    # e = field2[4]
    # f = field2[5]
    # g = field2[6]
    # h = field2[7]
fields3 = [[] for i in range(len(fields))]
fields4 = [[] for i in range(len(fields))]
fields5 = []
# #-----------------------------------------------------------------------


for j in  range (0, len(fields2), 1):
    for i in range (1, len(fields), 1):
        if float(fields2[j][0]) <= float(fields[i][3]) and float(fields2[j][0]) >= float(fields[i][3]) - 300:
            fields3[i].append(fields2[j][0])

print(fields3)

# ------------------------------for internal----------------------------------------------

for j in  range (0, len(fields2), 1):
    for i in range (1, len(fields), 1):
        if float(fields2[j][0]) > float(fields[i][3]) and float(fields2[j][0]) <= float(fields[i][4]):
            fields4[i].append("internal" + "_" + fields2[j][0])



# --------------------------------------------------------------------------------------------
# for i in range(len(fields)):
#     fields5.append(fields[i] + fields3[i] + fields4[i])
#
# for line in fields5:
#     print(*line)

for i in range(len(fields)):
    print(*fields[i] +  fields3[i])


# for i in range(len(fields)):
#     print(*fields[i] +  fields4[i])












#--------------------for negative strand------------------------------------------
#---------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------
F = open("operon_list_minus.txt")
f = F.readlines()
fields=[]

for line in f:
    field = line.strip("\n").split("\t")
    fields.append(field)
    a = field[0]
    b = field[1]
    c = field[2]
    d = field[3]
    e = field[4]
    f = field[5]
    g = field[6]
    # h = field[7]


F2 = open("_6_TSS_minus")
f2 = F2.readlines()
fields2=[]
for line in f2:
    field2 = line.strip("\n").split()
    fields2.append(field2)
    a = field2[0]
    b = field2[1]
    # c = field2[2]
    # d = field2[3]
    # e = field2[4]
    # f = field2[5]
    # g = field2[6]
    # h = field2[7]

fields3 = [[] for i in range(len(fields))]
fields4 = [[] for i in range(len(fields))]
fields5 = []


#-----------------------------------------------------------------------------------
for j in  range (0, len(fields2), 1):
    for i in range (1, len(fields), 1):
        if int(fields2[j][0]) >= int(fields[i][4]) and int(fields2[j][0]) <= int(fields[i][4]) + 300: # does not matter if we use float or int in this case
            fields3[i].append(fields2[j][0])




# ------------------------------for internal----------------------------------------------

# for j in  range (0, len(fields2), 1):
#     for i in range (1, len(fields), 1):
#         if float(fields2[j][0]) >= float(fields[i][3]) and float(fields2[j][0]) < float(fields[i][4]):
#             fields4[i].append("internal" + "_" + fields2[j][0])



#--------------------------------------------------------------------------------------------
# for i in range(len(fields)):
#     fields5.append(fields[i] + fields3[i] + fields4[i])

# for line in fields5:
#     print(*line)

for i in range(len(fields)):
    print(*fields[i] + fields3[i])
    # print(*fields[i] + fields3[i] , len(fields3[i]))

#
# for i in range(len(fields)):
#     print(*fields[i] +  fields4[i])
    # print(*fields[i] +  fields4[i] + len(fields4[i]))
