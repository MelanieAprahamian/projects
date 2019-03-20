#This script writes out the x and y coordinates for plotting the ROC curve for the known binders and decoy structures for troponin C.
#This script also spits out the area under the curve and EF.

import numpy as np
import matplotlib.pyplot as plt

########################################################
num_actives = 7 #number active ligands
num_decoys = 1000 #number decoy ligands
cutoff = 40 #cutoff for EF calculation

active_compounds = ['bepridil','w7','dfbp-o','trifluoperazine','pimobendan','nsc147866','levosimendan']

score = []
name = []

#determine how many poses were reported
#decoys:
with open('decoys.rept') as decoy_file:
    for line in decoy_file:
        if "decoy" in line and "/nfs" not in line:
            name.append(line.split()[1])
            score.append(line.split()[3])
            
#active ligands:
with open('actives.rept') as active_file:
    for line in active_file:
        if any(x in line for x in active_compounds):
            name.append(line.split()[1].split('_')[0])
            score.append(line.split()[3])

#combine name and score lists and order them into a single list, 'sorted_list'
zipped = zip(score,name)
sorted_list = sorted(zipped)

#split the tuples in the sorted list
score, name = zip(*sorted_list)

#ensure that if multiple structures have the same score, the active is put first
ordered_score = []
ordered_names = []

for i in range(len(score)):
    if i == 1:
        ordered_score.append(score[i])
        ordered_names.append(name[i])

    if i > 0:
        if score[i] != score[i-1]:
            position = i
            score_change = position
        if score_change == position and name[i].split('_', 1)[0] in active_compounds:
            ordered_score.insert(position, score[i])
            ordered_names.insert(position, name[i])
        if name[i].split('_', 1)[0] not in active_compounds:
            ordered_score.append(score[i])
            ordered_names.append(name[i])

#go through list and only keep top scoring of each active ligand and decoy
decoys = []

for i in range(1,num_decoys+1):
    decoys.append('decoy' + str(i))

test_name = []
new_name = []
new_score = []

for i in range(len(ordered_names)):
    if ordered_names[i] not in test_name:
        new_name.append(ordered_names[i])
        new_score.append(ordered_score[i])
        if any(activename in ordered_names[i] for activename in active_compounds) or any(decoyname in ordered_names[i] for decoyname in decoys):
            test_name.append(ordered_names[i])            

#add entries to the end of the lists to account for the decoys or actives that did not get scored
#the total length of both lists needs to be 1007
count_actives = 0
count_decoys = 0

total = num_actives + num_decoys

for i in test_name:
    if i in active_compounds:
        count_actives += 1
    else:
        count_decoys += 1

if len(new_score) < total:
    if count_actives < num_actives:
        for k in range(total - len(new_score)-(num_actives - count_actives)):
            new_name.append('decoy')
            new_score.append('0')
        for v in range(num_actives - count_actives):
            new_name.append('active')
            new_score.append('0')
    else:
        for k in range(total - len(new_score)):
            new_name.append('decoy')
            new_score.append('0')                
    
#ROC CURVE
#iterate through the sorted list
#start with point (0,0)
#if element in z is an active ligand (known binder) add 1 to y, else add 1 to x
count_x = 0
count_y = 0
x = []
y = []
area = 0.0
endpoint = 0.0
height = 0.0
m = 0
active_count = 0

#iterate through the binders/decoys and adjust the x,y coordinates accordingly and calculate area under curve
for n in range(len(new_score)):
    #add to y coordinate if active
    if new_name[n] in active_compounds:
        count_y += 1
        y.append(float(count_y)/float(num_actives))
        x.append(float(count_x)/float(num_decoys))
        active_count += 1
        if active_count == 1:
            length = float(x[n])
        if active_count > 1:
            length = float(x[n]) - float(x[m])
        height = float(y[m])
        area += (float(length) * float(height))
        m = n
    #add to x coordinate if non-active
    else:
        count_x += 1
        x.append(float(count_x)/float(num_decoys))
        y.append(float(count_y)/float(num_actives))

length = float(1.0) - float(x[m])
height = float(y[m])
area += (float(length) * float(height))

plt.figure()
plt.plot(x,y,linestyle='-',markersize=0,label='AUC = ' + str(area))
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.ylim(0,1)
plt.xlim(0,1)
plt.legend(loc='lower right')
plt.show()

#EF calculation
na = 0
for n in range(0,cutoff):
    if new_name[n] in active_compounds:
        na += 1
ef = (float(na)/float(cutoff))/(float(num_actives)/(float(num_actives)+float(num_decoys)))
