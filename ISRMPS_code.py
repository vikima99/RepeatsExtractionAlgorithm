import timeit
import re

#################################### Collecting all sequences ################################

def collect_All_sequence(path):
    file=open(path,'r')
    f=file.readlines()
    seq=[]
    s=''
    for line in f:       
        if line[0]=='>':
            seq.append(s)
            s=''
        elif line[0]!='>':
            s+=line.strip('\n')
    seq.append(s)
    seq=seq[1:]    
    return seq  

######################################## SEQUENCE and SEQUENCE_ID FINDER ########################################

def collect_All_ID(path):
    file=open(path,'r')
    f=file.readlines()
    dict_id={}
    count=0
    for line in f:       
        if line[0]=='>':
            count+=1
            regex=re.search(r"\|(|.+|)\|",line)
            match=regex.group(1)
            dict_id[count]=match    
        elif line[0]!='>':
            pass
            # s+=line.strip('\n')
    return dict_id

################-------------------------------------------------------------------------------------################

def sequence_finder(unique_sequence,index):
    spl_char_count=0
    start=0
    end=len(unique_sequence[0])
    if index>=start and index <=end:
        spl_char_count+=1
        correct_index=index-start
        return 1,correct_index+1
    for i in range(1,len(unique_sequence)):
        spl_char_count+=1
        if spl_char_count<=len(spl_char):
            start=end+1
            end=end+len(unique_sequence[i])+1
            if index>=start and index <=end:
                correct_index=index-start
                return i+1,correct_index+1           
        else:
            start=end+3
            end=end+len(unique_sequence[i])+3
            if index>=start and index <=end:
                correct_index=index-start
                return i+1,correct_index-1

#######################################  Initial Hash Computation ########################################

def hash(Amino_acids,spl_char,max_num):
    num=max_num+3
    dict={}
    for i in range(len(Amino_acids)):
        dict[Amino_acids[i]]=ord(Amino_acids[i])
    
    for i in range(len(spl_char)):
        dict[spl_char[i]]=num
        num=num+1
    return dict
# coded=hash(Amino_acids,spl_char,max_num)


################# Combining All sequences into one using unique spl characters ##############################

def Join_string(unique_sequence,spl_char):
    if len(unique_sequence)<=len(spl_char):
        vik=''
        for i in range(len(unique_sequence)):
            vik+=unique_sequence[i]
            vik+=spl_char[i]
        vik=vik[:-1]
    else:
        vik=''
        nume_=0
        for i in range(len(unique_sequence)):
            vik+=unique_sequence[i]
            ## apply condition ##
            if (i+1)<len(spl_char):
               vik+=spl_char[i]
            else:
                nume_+=1
                join_char=Kth_combination(spl_char,nume_)
                vik+=join_char


        # vik=vik[:-1]

    return vik


######################################################################################################################

def rabin_karp_2(step,nums,l,S):
    p=2**63 - 1
    cur_hash=0
    for i in range(step):
        cur_hash = (cur_hash * l + nums[i]) % p
    hashes = [cur_hash]    
    max_pow = pow(l,step, p)
    for i in range(step, len(S)):
        cur_hash = (l*cur_hash-nums[i-step]*max_pow + nums[i]) % p
        hashes.append(cur_hash)

    dict2={}
    dict2_={}
    
    for i in range(len(hashes)):
        
        if dict2.get(hashes[i],'nf')=='nf':
            dict2[hashes[i]]=1
            dict2_[hashes[i]]=[i]
        else:
            dict2[hashes[i]]+=1
            vik=dict2_[hashes[i]]
            vik.append(i)
            dict2_[hashes[i]]=vik

    dict3={}
    for key in dict2:
        if dict2[key]>1:
            dict3.update({key:[dict2[key],dict2_[key],step]})
    return dict3
   
#################### ###################################### Rabin_Karp_Algorithm ##################################################################

def Rabin_Karp(mid_value,nums,p,l,S):                    # l is no of thing possible
    cur_hash=0
    for i in range(mid_value):
        cur_hash = (cur_hash * l + nums[i]) % p
    hashes = {cur_hash}
    pos = -1
    max_pow = pow(l, mid_value, p)
    for i in range(mid_value, len(S)):
        cur_hash = (l*cur_hash-nums[i-mid_value]*max_pow + nums[i]) % p
        if cur_hash in hashes:
            pos = i + 1 - mid_value
        hashes.add(cur_hash)
    return pos

########################################### Rabin-Karp Pattern searching Algorithm ##################################################

def getZarr(string, z):
    n = len(string)
    l, r, k = 0, 0, 0
    for i in range(1, n):
        if i > r:
            l, r = i, i
            while r < n and string[r - l] == string[r]:
                r += 1
            z[i] = r - l
            r -= 1
        else:
            k = i - l
            if z[k] < r - i + 1:
                z[i] = z[k]
            else:
                l = i
                while r < n and string[r - l] == string[r]:
                    r += 1
                z[i] = r - l
                r -= 1

#-------------------------------------------------------------------------------------------------------------#

def search(text, pattern):
    concat = pattern + spl_char[-1] + text
    l = len(concat)
    z = [0] * l
    getZarr(concat, z)
    list=[]
    for i in range(l):
        if z[i] == len(pattern):
            # # print("Pattern found at index",
            #           i - len(pattern) - 1)
            num_=i-len(pattern)-1
            list.append(num_)
    return list
########################################## searchin regex #############################################

def search_new(text,pattern):
    posi=[]
    for m in re.finditer(pattern, text):
        posi.append(m.start())
    return posi

####################################### Kth combination ##################################################

def Kth_combination(iterable,index):
    pool = tuple(iterable)
    n = len(pool)
    r=3
    if r < 0 or r > n:
        raise ValueError
    c = 1
    k = min(r, n-r)

    for i in range(1, k+1):
        c = c * (n - k + i) // i

    if index < 0:
        index += c

    if index < 0 or index >= c:
        raise IndexError

    result = []
    while r:
        c, n, r = c*r//n, n-1, r-1
        while index >= c:
            index -= c
            c, n = c*(n-r)//n, n-1
        result.append(pool[-1-n])

    vik=''
    for char in result:
        vik+=char
    return vik


###################################### SEARCHING LONGEST DUPLICATE SUBSTRING ##################################

def LongestDupliString(S,l):
    p=2**63 - 1
    
    low, high = 0, len(S)-1
    
    start = 0
    coded=hash(Amino_acids,spl_char,max_num)
    # print(S)
    # for c in S:
    #   print(f"{c} {coded[c]-50}")

    nums = [coded[c]-50 for c in S]
    
    while low <= high:
        mid_value = (low+high) // 2
        pos = Rabin_Karp(mid_value,nums,p,l,S)
        if pos == -1:                   # no matching strings found
            high = mid_value - 1
        else:
            start = pos
            low = mid_value+ 1
    # print(f"longest duplicate string is : {S[start:start+low-1]}")
    
    max_len=len(S[start:start+low-1])
    # print(f"Its length : {max_len}")
    
    loc_=rabin_karp_2(max_len, nums, l, S)

    return loc_,max_len

###############################################*** RECURSIVELY PATTERN FINDING **################################################

def New_sentence_creation(seq,pattern_length,index,spl_char,remove_index):

    if index<=(len(spl_char)-1):   
        i=0
        factor=pattern_length-1
        for value in remove_index:
            if index<=(len(spl_char)-1):
                value=value-factor*i
                seq=seq[:value]+spl_char[index]+seq[(value+pattern_length):]
                # print(seq)
                i+=1
                index+=1
            else:
                break
    else:
        i=0
        factor=pattern_length-3
        for value in remove_index:
            value=value-factor*i
            if len(unique_sequence)<=len(spl_char):
                combi_index=index-(len(spl_char)-1)
            else:
                combi_index=index-(len(unique_sequence)-len(spl_char)-5)
            join_char=Kth_combination(spl_char,combi_index)
            seq=seq[:value]+join_char+seq[(value+pattern_length):]
            i+=1
            index+=1
    return seq,index+1

#############################################------------------------------#####################################################

def Recursion_support(max_occur,max_val,loc_,spl_char,spl_char_index,seq,Longest_pattern_length,pat,max_):
    if Longest_pattern_length>display_threshold:
        no_of_pattern=list(loc_.values())
        for p in no_of_pattern:
            length_of_pattern=p[2]
            first_index=p[1][0]
            pat_sign=seq[first_index:(first_index+length_of_pattern)]
            if Amino_acids.count(pat_sign[-1])==0:
                pat_sign=pat_sign[:-2]
            if len(pat_sign)>display_threshold:
                if pat_sign not in pat:
                    a,b=flushing_pattern_found(pat_sign,S,file_seq,max_occur)
                    if a!='Discard it!':
                        max_occur=b
                        if max_occur>max_val:
                            max_val=max_occur
                            max_.add(max_val)
                        pat.add(pat_sign)
                        file_seq.write(f"Sl.No. of the repeat: {len(list(pat))+1}")
                        file_seq.write('\n')
        remove_index=[]
        for p in no_of_pattern:
            posi=p[1]
            for k in range(len(posi)-1):
                # print(posi[k])
                remove_index.append(posi[k])
        pattern_length=length_of_pattern
        index=spl_char_index
        seq,spl_char_index=New_sentence_creation(seq,pattern_length,index,spl_char,remove_index) 
        loc_,Longest_pattern_length=LongestDupliString(seq,l)
        
    return max_occur,max_val,loc_,spl_char,spl_char_index,seq,Longest_pattern_length,pat,max_




#############################################-------------------------------####################################################

def string_creation_recursion(max_occur,max_val,loc_,spl_char,spl_char_index,seq,Longest_pattern_length,pat=set(),max_=set()):
    while (Longest_pattern_length>display_threshold):
        no_of_pattern=list(loc_.values())
        for p in no_of_pattern:
            length_of_pattern=p[2]
            first_index=p[1][0]
            pat_sign=seq[first_index:(first_index+length_of_pattern)]
            if Amino_acids.count(pat_sign[-1])==0:
                pat_sign=pat_sign[:-2]
            if len(pat_sign)>display_threshold:
                if pat_sign not in pat: 
                    a,b=flushing_pattern_found(pat_sign,S,file_seq,max_occur)
                    if a!='Discard it!':
                        max_occur=b
                        if max_occur>max_val:
                            max_val=max_occur
                            max_.add(max_val)
                        pat.add(pat_sign)
                        # print(f"Repeat : {pat_sign}")
                        # print('\n')
                        file_seq.write(f"Sl.No. of the repeat: {len(list(pat))+1}")
                        file_seq.write('\n')
        remove_index=[]
        for p in no_of_pattern:
            posi=p[1]
            for k in range(len(posi)-1):
                remove_index.append(posi[k])
        pattern_length=length_of_pattern
        index=spl_char_index
        seq,spl_char_index=New_sentence_creation(seq,pattern_length,index,spl_char,remove_index) 
        loc_,Longest_pattern_length=LongestDupliString(seq,l)
        max_occur,max_val,loc_,spl_char,spl_char_index,seq,Longest_pattern_length,pat,max_=Recursion_support(max_occur,max_val,loc_,spl_char,spl_char_index,seq,Longest_pattern_length,pat,max_)
        # string_creation_recursion(max_occur,max_val,loc_,spl_char,spl_char_index,seq,Longest_pattern_length,pat=pat,max_=max_)
    return 'Completed Successfully !!',len(list(pat)),list(pat),max_

##########################################################  Displaying Code ############################################################

def flushing_pattern_found(pattern,S,file_seq,max_occur):
    text = S
    # posi=search(text,pattern)
    posi=search_new(text,pattern)
    if len(posi)>max_occur:
        max_occur=len(posi)
    if len(posi)==1 or len(posi)==0:
      return 'Discard it!','hello'
    else:
        file_seq.write(f"No. of amino acids in the repeat: {len(pattern)}")
        # file_seq.write(f"REPEAT > :Length : [{len(pattern)}]")
        file_seq.write('\n')
        file_seq.write(f"No. of occurrences : {len(posi)}")
        # file_seq.write('\n')
        file_seq.write('\n\n')
        file_seq.write('Repeat:')
        file_seq.write('\n')
        file_seq.write(pattern)
        file_seq.write('\n\n')
        # file_seq.write('\t')

        # print(f"REPEAT > :Length : {len(pattern)} :> {pattern}  : Ocurrences : {len(posi)}")
        for k in range(len(posi)):
            v=sequence_finder(unique_sequence,posi[k])
            v=list(v)
            file_seq.write(f"Sequence ID : {id_hash[v[0]]} | No. of amino acids in this sequence: [{len(unique_sequence[v[0]-1])}] | Position : {[(v[1]),(v[1]+len(pattern)-1)]} |")
            file_seq.write('\n')
            # print(f"Query No > {v[0]} and Position : {[(v[1]-1),(v[1]-1+len(pattern)-1)]}")
        file_seq.write('\n')
        file_seq.write('\n')
        # print('\n')
        file_seq.write('__________________________________________________________________________________________________________________________')
        file_seq.write('\n\n')
        return 'Real pattern',max_occur
   
    
    


################################################# ** Main code **  ######################################################

######## -------path to spl char -----------###
# path_spl=r'C:\\web_server\\special_char_3500.txt'

path_spl=r"/mnt/c/web_server/special_char_3500.txt"
sp=open(path_spl)
spl=sp.readlines()
spl_char=[]
for j in spl:
    get=j.strip('\n')
    spl_char.append(get)
spl_char=spl_char[:2502]


Amino_acids='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
max_num=91


#######----------path to sequence file ------------------###

# import os
# path_initial=r"/mnt/c/web_server/Sequence_Files/FINALORGANS"
# path_write_i=r"/mnt/c/web_server/Repeats_Files/All_Organ_Output"
# path_write_i=r"/mnt/c/web_server/Repeats_Files/TB_repeats/GLYCOSYLATIONBLASTFORMATCASESTUDY.fasta"

path_write_i=r"/mnt/c/web_server/Test_input_seq_repeats.fasta"


# path_write_i=r"c\\web_server\\sample_today.fasta"


######################################--------------------- ###########################################################################

# display_threshold=input("Display Threshold of the repeat ")
# display_threshold=int(display_threshold)

display_threshold=4

# path_ini=r"/mnt/c/web_server/Genome/chromosome_files"
# path_ini=r"/mnt/c/web_server/Sequence_Files/GLYCOSYLATIONBLASTFORMATCASESTUDY.txt"
path_ini=r"/mnt/c/web_server/Test_case.fasta"
# path_ini=r"c\\web_server\\sample_anirban.txt"

# path_ini=r'/mnt/c/Users/vikim/OneDrive/Desktop/Monday_saperate_work/prtn_sequence/Covid_case_study.fasta'

# path_ini=r"/mnt/c/web_server/Two_example_seq.fasta"
path_to_write=path_write_i
# print(file)
final_path=path_ini
max_occur=0
max_val=max_occur
start = timeit.default_timer()
print("file is getting ready!!")
# print('\n')
id_hash=collect_All_ID(final_path)
unique_sequence=collect_All_sequence(final_path)
l=26+len(unique_sequence)+len(spl_char)
S=Join_string(unique_sequence,spl_char)
v,Longest_pattern_length=LongestDupliString(S,l)
seq=S
spl_char_index=len(unique_sequence)-1
loc_=v

##########################################---------------------------------main writing code ---------------------------------#################

file_seq=open(path_to_write,'w')
file_seq.write(f'ISRMPS: A method to locate Identical Sequence Repeats across Multiple Protein Sequences.')
file_seq.write('\n\n')
file_seq.write('________________________________________________________________________________________________________________________')
file_seq.write('\n\n')
file_seq.write(f"Sl.No. of repeat : {1}")
file_seq.write('\n')
Message,count,min_length,max_=string_creation_recursion(max_occur,max_val,loc_,spl_char,spl_char_index,seq,Longest_pattern_length,pat=set(),max_=set())
mini_=[]
for seq_ in min_length:
    mini_.append(len(seq_))
min_pat=min(mini_)
file_seq.close()
end = timeit.default_timer()
with open(path_to_write) as f1:
    lines = f1.readlines()
N=1
with open(path_to_write, 'w') as file_seq:
    file_seq.writelines(lines[:-N])
    file_seq.write("Summary of the output:")
    file_seq.write('\n')
    file_seq.write('__________________________________________________________________________________________________________________________________')
    file_seq.write('\n')
    file_seq.write('\n')
    file_seq.write(f"Total No. of input sequences : {len(unique_sequence)}")
    file_seq.write('\n')
    file_seq.write(f"Total No. of repeats found : {count}")
    file_seq.write('\n')
    file_seq.write('\n')
    file_seq.write(f"Minimum No. of amino acids in repeat: {min_pat}")
    file_seq.write('\n')
    file_seq.write(f"Maximum No. of amino acids in repeat: {Longest_pattern_length}")
    file_seq.write('\n')
    file_seq.write('\n')
    file_seq.write(f'Maximum No. of ocurrences: {max(max_)}')
    file_seq.write('\n')
    file_seq.write('Minimum No. of ocurrences: 2')
    file_seq.write('\n')
    file_seq.write('\n')
    file_seq.write(f"Time taken : {end-start} seconds OR : {(end-start)/60} minutes")
    print(f"Time taken : {end-start} seconds OR : {(end-start)/60} minutes")
    file_seq.write('\n')
    file_seq.write('\n')
    file_seq.close()


    
########################---------------------################################







