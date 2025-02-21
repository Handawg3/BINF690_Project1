#!/usr/bin/env python

#import re, although I don't think I ever used it
import re
#load in the exclude list
F = open('c:/Users/andre/OneDrive/Documents/BINF690/Project 1/exclude_list.tab')

#initialize a list variable and store the scaffold names of the exclude list in it
holder = []
for line in F:
    line = line.strip()
    holder.append(line.split('\t')[0])

#open the output file
output = open('c:/Users/andre/OneDrive/Documents/BINF690/Project 1/afterscaffold.fasta','w')

#load in the fasta file 
S = open('c:/Users/andre/OneDrive/Documents/BINF690/Project 1/Scaffolds.fasta')

#load in each line of S strip it, then if it starts with a >, so a scaffold line,
#and the scaffold name in the scaffold file is not in the holder list then output the line
#and the next line which is the sequence.
for line in S:
    line = line.strip()
    if line.startswith('>'):
        if line[1:] not in holder:
            output.write(line+'\n')
            output.write(S.readline())

#close files used
output.close()
F.close()
S.close()

#open trim list
T = open('c:/Users/andre/OneDrive/Documents/BINF690/Project 1/trim_list.tab')

#separate each line by the tabs and store each variable in the line into a list
name = []
length = []
span = []
info = []
for line in T:
    if line.startswith('s'):
        line = line.strip()
        name.append(line.split('\t')[0])
        length.append(line.split('\t')[1])
        span.append(line.split('\t')[2])
        info.append(line.split('\t')[3])

#open an output file for question 2 and open the scaffold file after exclusion
P = open('c:/Users/andre/OneDrive/Documents/BINF690/Project 1/afterscaffold.fasta')
output = open('c:/Users/andre/OneDrive/Documents/BINF690/Project 1/aftertrim_Q2.fasta','w')        

#establish a threshold for minimum sequence length
threshold = 200
for line in P:
    #if line starts with >, it is the scaffold, then strip it of \n for easier matching
    if line.startswith('>'):
        line = line.strip()
        #if the line in the scaffold file is in the name list of the trim list
        if line[1:] in name:
            #store the value of the index of the match in a variable since the trim list is not in order
            idx = name.index(line[1:])
            seqline = P.readline().strip()
            #if the sequence at idx we are working with doesn't have the mito... info, we trim it
            if info[idx] != 'mitochondrion-not_cleaned':
                #split the span of the current sequence at the comma, then store in a tuple after then
                #splitting at the "..", and establish a counter to work through the list within the list
                #from the first split, we can do this since lists are muttable
                span[idx] = span[idx].split(',')
                counter = 0
                for i in span[idx]:
                    span[idx][counter] = tuple(span[idx][counter].split('..'))
                    counter += 1
                #set a boolean value to true for ease of working through loop, set a count to zero and initialize trim
                x = True
                trim = ''
                count = 0
                #go through each span set of the corresponding sequence
                for t in span[idx]:
                    #if the first span starts with 1, we change the trim to the new sequence without the first span,
                    #and set boolean to false so the next time through we end in the middle block
                    if t[0] == '1':
                        trim = seqline[(int(t[1])-1):]
                        last = int(t[1])
                        x = False
                    #this block will only go into if one of the others has already been used first because the boolean should be true until that happens
                    elif x == False:
                        #trim will equal from the second value in the last tuple to the first value in the new one minus 1 for the indexing
                        trim = seqline[last:((int(t[0]))-1)]
                        #only output the sequence if its length is higher than the threshold and it has been stripped of N's on the ends
                        if len(trim) > threshold:
                            output.write(line+f'_{count}\n')
                            trim = trim.strip('N')  
                            output.write(trim+'\n')
                        last = int(t[1])
                    #this block will only hit if the first value of the first tuple is not 1, to make sure from the beginning to the first value in the 
                    #tuple is printed
                    else:
                        #have to start count at 1 here or it will print the subsequence with a _0
                        count = 1
                        trim = seqline[:(int(t[0])-1)] 
                        #only output the sequence if its length is higher than the threshold and it has been stripped of N's on the ends
                        if len(trim) > threshold:
                            output.write(line+f'_{count}\n')
                            trim = trim.strip('N')  
                            output.write(trim+'\n')
                        x = False
                        last = int(t[1])
                    #increase count for the naming of the subsequences
                    count += 1
                #else loop on the for block to check if the second value of the last tuple is the same as the length of the sequence
                #if so then that is trimmed and nothing needs to happen
                else:
                    if last == len(seqline):
                        last = 0
                    #if it isn't the same from the last value, to the end is the last subsequence
                    else:
                        trim = seqline[last:]
                        #only output the sequence if its length is higher than the threshold and it has been stripped of N's on the ends
                        if len(trim) > threshold:
                            output.write(line+f'_{count}\n')
                            trim = trim.strip('N')
                            output.write(trim+'\n')
            #this is the else to the mito.. if condition, if it is mito.., print as normal and 
            #only output the sequence if its length is higher than the threshold and it has been stripped of N's on the ends
            else:
                if len(seqline) > threshold:
                    output.write(line+'\n')
                    output.write(seqline+'\n')
        #The else if the scaffold isn't in the name list, print as normal and 
        #only output the sequence if its length is higher than the threshold and it has been stripped of N's on the ends
        else:
            seqline = P.readline()
            if len(seqline) > threshold:
                output.write(line+'\n')
                output.write(seqline)

#close files used
P.close()
output.close()
T.close()

#open all the files needed for part 3
T = open('c:/Users/andre/OneDrive/Documents/BINF690/Project 1/trim_list.tab')
P = open('c:/Users/andre/OneDrive/Documents/BINF690/Project 1/afterscaffold.fasta')
output = open('c:/Users/andre/OneDrive/Documents/BINF690/Project 1/aftertrim_Q3.fasta','w')

#separate each line by the tabs and store each variable in the line into a list
name = []
length = []
span = []
info = []
for line in T:
    if line.startswith('s'):
        line = line.strip()
        name.append(line.split('\t')[0])
        length.append(line.split('\t')[1])
        span.append(line.split('\t')[2])
        info.append(line.split('\t')[3])

#establish a threshold for minimum sequence length
threshold = 200
for line in P:
    #if line starts with >, it is the scaffold, then strip it of \n for easier matching
    if line.startswith('>'):
        line = line.strip()
        #if the line in the scaffold file is in the name list of the trim list
        if line[1:] in name:
            #store the value of the index of the match in a variable since the trim list is not in order
            idx = name.index(line[1:])
            seqline = P.readline().strip()
            #split the span of the current sequence at the comma, then store in a tuple after then
            #splitting at the "..", and establish a counter to work through the list within the list
            #from the first split, we can do this since lists are muttable
            span[idx] = span[idx].split(',')
            counter = 0
            for i in span[idx]:
                span[idx][counter] = tuple(span[idx][counter].split('..'))
                counter += 1
            #set a boolean value to true for ease of working through loop, set a count to zero and initialize trim
            x = True
            trim = ''
            count = 0
            #go through each span set of the corresponding sequence
            for t in span[idx]:
                #if the first span starts with 1, we change the trim to the new sequence without the first span,
                #and set boolean to false so the next time through we end in the middle block
                if t[0] == '1':
                    trim = seqline[(int(t[1])-1):]
                    last = int(t[1])
                    x = False
                #this block will only go into if one of the others has already been used first because the boolean should be true until that happens
                elif x == False:
                    #trim will equal from the second value in the last tuple to the first value in the new one minus 1 for the indexing
                    trim = seqline[last:((int(t[0]))-1)]
                    #only output the sequence if its length is higher than the threshold and it has been stripped of N's on the ends
                    if len(trim) > threshold:
                        output.write(line+f'_{count}\n')
                        trim = trim.strip('N')  
                        output.write(trim+'\n')
                    last = int(t[1])
                #this block will only hit if the first value of the first tuple is not 1, to make sure from the beginning to the first value in the 
                #tuple is printed
                else:
                    #have to start count at 1 here or it will print the subsequence with a _0
                    count = 1
                    trim = seqline[:(int(t[0])-1)] 
                    #only output the sequence if its length is higher than the threshold and it has been stripped of N's on the ends
                    if len(trim) > threshold:
                        output.write(line+f'_{count}\n')
                        trim = trim.strip('N')  
                        output.write(trim+'\n')
                    x = False
                    last = int(t[1])
                #increase count for the naming of the subsequences
                count += 1
            #else loop on the for block to check if the second value of the last tuple is the same as the length of the sequence
            #if so then that is trimmed and nothing needs to happen
            else:
                if last == len(seqline):
                    last = 0
                #if it isn't the same from the last value, to the end is the last subsequence
                else:
                    trim = seqline[last:]
                    #only output the sequence if its length is higher than the threshold and it has been stripped of N's on the ends
                    if len(trim) > threshold:
                        output.write(line+f'_{count}\n')
                        trim = trim.strip('N')
                        output.write(trim+'\n')
        #The else if the scaffold isn't in the name list, print as normal and 
        #only output the sequence if its length is higher than the threshold and it has been stripped of N's on the ends
        else:
            seqline = P.readline()
            if len(seqline) > threshold:
                output.write(line+'\n')
                output.write(seqline)

#close files used
P.close()
output.close()
T.close()