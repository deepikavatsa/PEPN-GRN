
## compute Score, p_AUROC,AUROC, p_AUPR, AUPR using dream4 scoring method

from dreamtools import D4C2
s = D4C2()

#######################################  KMEANS	####################################

f = open('dream4_scores_kmeans.txt','w')

## computation for 10-node Dream4 data sets
f.write("\n==============================================================================\n")
f.write("\n\t\tPEPN-GRN-v3 (kmeans Discretized Dream4 10-node data sets)\n")
f.write("\n==============================================================================\n")

filename1="all_edges_10_1_kmeans.txt"
filename2="all_edges_10_2_kmeans.txt"
filename3="all_edges_10_3_kmeans.txt"
filename4="all_edges_10_4_kmeans.txt"
filename5="all_edges_10_5_kmeans.txt"
scores = s.score(filename1,"10")
f.write("\n10_1: \n")
f.write(str(scores))

scores = s.score(filename2,"10")
f.write("\n10_2: \n")
f.write(str(scores))

scores = s.score(filename3,"10")
f.write("\n10_3: \n")
f.write(str(scores))

scores = s.score(filename4,"10")
f.write("\n10_4: \n")
f.write(str(scores))

scores = s.score(filename5,"10")
f.write("\n10_5: \n")
f.write(str(scores))

scores = s.score([filename1,filename2,filename3,filename4,filename5],"10")
f.write("\n\nFor all 5 networks: \n")
f.write(str(scores))


## computation for 100-node Dream4 data sets
f.write("\n\n==============================================================================\n")
f.write("\n\t\tPEPN-GRN (kmeans 3-bin Discretized Dream4 100-node data sets)\n")
f.write("\n==============================================================================\n")

filename1="all_edges_100_1_kmeans.txt"
filename2="all_edges_100_2_kmeans.txt"
filename3="all_edges_100_3_kmeans.txt"
filename4="all_edges_100_4_kmeans.txt"
filename5="all_edges_100_5_kmeans.txt"
scores = s.score(filename1,"100")
f.write("\n100_1: \n")
f.write(str(scores))

scores = s.score(filename2,"100")
f.write("\n100_2: \n")
f.write(str(scores))

scores = s.score(filename3,"100")
f.write("\n100_3: \n")
f.write(str(scores))

scores = s.score(filename4,"100")
f.write("\n100_4: \n")
f.write(str(scores))

scores = s.score(filename5,"100")
f.write("\n100_5: \n")
f.write(str(scores))

scores = s.score([filename1,filename2,filename3,filename4,filename5],"100")
f.write("\n\nFor all 5 networks: \n")
f.write(str(scores))

f.close()

########################################	EFD	####################################
f = open('dream4_scores_efd.txt','w')

## computation for 10-node Dream4 data sets
f.write("\n==============================================================================\n")
f.write("\n\t\tPEPN-GRN (efd Discretized Dream4 10-node data sets)\n")
f.write("\n==============================================================================\n")

filename1="all_edges_10_1_efd.txt"
filename2="all_edges_10_2_efd.txt"
filename3="all_edges_10_3_efd.txt"
filename4="all_edges_10_4_efd.txt"
filename5="all_edges_10_5_efd.txt"
scores = s.score(filename1,"10")
f.write("\n10_1: \n")
f.write(str(scores))

scores = s.score(filename2,"10")
f.write("\n10_2: \n")
f.write(str(scores))

scores = s.score(filename3,"10")
f.write("\n10_3: \n")
f.write(str(scores))

scores = s.score(filename4,"10")
f.write("\n10_4: \n")
f.write(str(scores))

scores = s.score(filename5,"10")
f.write("\n10_5: \n")
f.write(str(scores))

scores = s.score([filename1,filename2,filename3,filename4,filename5],"10")
f.write("\n\nFor all 5 networks: \n")
f.write(str(scores))


## computation for 100-node Dream4 data sets
f.write("\n\n==============================================================================\n")
f.write("\n\t\tPEPN-GRN (efd 3-bin Discretized Dream4 100-node data sets)\n")
f.write("\n==============================================================================\n")

filename1="all_edges_100_1_efd.txt"
filename2="all_edges_100_2_efd.txt"
filename3="all_edges_100_3_efd.txt"
filename4="all_edges_100_4_efd.txt"
filename5="all_edges_100_5_efd.txt"
scores = s.score(filename1,"100")
f.write("\n100_1: \n")
f.write(str(scores))

scores = s.score(filename2,"100")
f.write("\n100_2: \n")
f.write(str(scores))

scores = s.score(filename3,"100")
f.write("\n100_3: \n")
f.write(str(scores))

scores = s.score(filename4,"100")
f.write("\n100_4: \n")
f.write(str(scores))

scores = s.score(filename5,"100")
f.write("\n100_5: \n")
f.write(str(scores))

scores = s.score([filename1,filename2,filename3,filename4,filename5],"100")
f.write("\n\nFor all 5 networks: \n")
f.write(str(scores))

f.close()

########################################	EWD	####################################

f = open('dream4_scores_ewd.txt','w')

## computation for 10-node Dream4 data sets
f.write("\n==============================================================================\n")
f.write("\n\t\tPEPN-GRN (ewd 3-bin Discretized Dream4 10-node data sets)\n")
f.write("\n==============================================================================\n")

filename1="all_edges_10_1_ewd.txt"
filename2="all_edges_10_2_ewd.txt"
filename3="all_edges_10_3_ewd.txt"
filename4="all_edges_10_4_ewd.txt"
filename5="all_edges_10_5_ewd.txt"
scores = s.score(filename1,"10")
f.write("\n10_1: \n")
f.write(str(scores))

scores = s.score(filename2,"10")
f.write("\n10_2: \n")
f.write(str(scores))

scores = s.score(filename3,"10")
f.write("\n10_3: \n")
f.write(str(scores))

scores = s.score(filename4,"10")
f.write("\n10_4: \n")
f.write(str(scores))

scores = s.score(filename5,"10")
f.write("\n10_5: \n")
f.write(str(scores))

scores = s.score([filename1,filename2,filename3,filename4,filename5],"10")
f.write("\n\nFor all 5 networks: \n")
f.write(str(scores))


## computation for 100-node Dream4 data sets
f.write("\n\n==============================================================================\n")
f.write("\n\t\tPEPN-GRN (ewd Discretized Dream4 100-node data sets)\n")
f.write("\n==============================================================================\n")

filename1="all_edges_100_1_ewd.txt"
filename2="all_edges_100_2_ewd.txt"
filename3="all_edges_100_3_ewd.txt"
filename4="all_edges_100_4_ewd.txt"
filename5="all_edges_100_5_ewd.txt"
scores = s.score(filename1,"100")
f.write("\n100_1: \n")
f.write(str(scores))

scores = s.score(filename2,"100")
f.write("\n100_2: \n")
f.write(str(scores))

scores = s.score(filename3,"100")
f.write("\n100_3: \n")
f.write(str(scores))

scores = s.score(filename4,"100")
f.write("\n100_4: \n")
f.write(str(scores))

scores = s.score(filename5,"100")
f.write("\n100_5: \n")
f.write(str(scores))

scores = s.score([filename1,filename2,filename3,filename4,filename5],"100")
f.write("\n\nFor all 5 networks: \n")
f.write(str(scores))


f.close()
