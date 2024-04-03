import numpy as np

def get_migrating_tracts(ts,pop,ind,L,pop_id=-1):
    """
    get_migrating_tracts_ind returns for a single individual ind, all the semgents which introgressed from population pop.
    :ts: the tskit Tree Sequence
    :pop: we recover segments with an ancestry from population named pop
    :ind: the indivual we are interested in recovering segments with population 'pop' ancestry
    :pop_id: an added possibility to give as an input the id of the population instead of its name
    :return: a list of segments with ancestry pop
    """ 
    if pop_id==-1:  #We recover the id of the pop we are interested in, except if it was directly given as an input
        pop_id = [p.id for p in ts.populations() if p.metadata['name']==pop][0]
    mig = ts.tables.migrations
    migrating_tracts = [] #The list of segments with ancestry pop.
    
    for tree in ts.trees(): #For each tree (which corresponds to a subpart of the simulated genome) we check if that segment has an ancestry in pop. If so we had it to our list of segments.
        u = ind 
        while u != tskit.NULL: #We loop through all the ancestors of ind in the tree and see if one of them migrated in pop. If so the current segment should be added to our list. 
            migs = np.where(mig.node == u)[0] #All the migration events of node u.
            for cur_mig in migs:
                cur_mig = mig[cur_mig]
                if(cur_mig.dest==pop_id and cur_mig.left<=tree.interval.left and cur_mig.right>=tree.interval.right): #We check if that migration is in 'pop' and if it is contained in the segment of our tree. If so we add the segment to our list.
                    flag=False
                    if(len(migrating_tracts)>0 and tree.interval.left==migrating_tracts[len(migrating_tracts)-1][1]): #A migration will likely overlap over multiple adjacent trees, we merge these adjacent segments if that is the ca
                        migrating_tracts[len(migrating_tracts)-1][1]=tree.interval.right                     
                    else:
                        migrating_tracts.append([tree.interval.left,tree.interval.right])
            u = tree.parent(u)
    migrating_tracts = clean_tracts(migrating_tracts,L)
    return migrating_tracts

#Return the first tract minus the second
def substract_tracts(tracts1,tracts2):
    maxi = 0
    res =[]
    for t in tracts1:
        if t[1]>maxi:
            maxi=t[1]
    for t in tracts2:
        if t[1]>maxi:
            maxi=t[1]
    b1=False
    b2=True
    inside=False
    start=-1
    for i in range(maxi+1):
        if (inTracts(i,tracts1)):
            b1=True
        else:
            b1=False
        if (inTracts(i,tracts2)):
            b2=False
        else:
            b2=True
        if (b1 and b2 and not inside):
            inside = True
            start = i
        elif( not(b1 and b2) and inside):
            inside = False
            res.append([start,i])

    if(b1 and b2 and inside):
        res.append([start,maxi])
    return res

#Check if the position is in the list of tracts
def inTracts(pos,tract):
    for i in range(len(tract)):
        if pos>=tract[i][0] and pos <= tract[i][1]:
            return True
    return False

#Create a list of tracts from HMM result
def get_HMM_tracts(seq):
    migrating_tracts = []
    maxi = 0
    for e in seq:
        if e>maxi:
            maxi=e
            
    for i in range(maxi+1):
        migrating_tracts.append([])
    start=0
    for i in range(1,len(seq)):
        if seq[i]!=seq[i-1]:
            migrating_tracts[seq[i-1]].append([start,i-1])
            start=i
    migrating_tracts[seq[len(seq)-1]].append([start,len(seq)-1])
    return migrating_tracts

# Input:
#    seq: The result of the HMM algorithm
#    tracts: The actual tracts, in the same order as the states. Eg. if the state 0 correspond to non Archaic, the first tracts should correspond to non Archaic tracts.
def confusionMatrix(seq, tracts):
    nbState =len(tracts)
    M = np.zeros((nbState,nbState),dtype=int)
    for j in range(nbState):
        for t in tracts[j]:
            for i in range(t[0],t[-1]):
                M[seq[i],j]+=1
    return M

def clean_tracts(tractInit,size):
    tract = np.copy(tractInit)
    tract = tract/size
    tract=tract.astype(int)
    flag = True
    while(flag):
        flag=False
        for i in range(len(tract)):
            for j in range(len(tract)):
                if not flag and tract[i,0]==tract[j,1]:
                    tract[j,1]=tract[i,1]
                    tract = np.delete(tract,i,0)
                    flag=True
    flag = True
    while(flag):
        flag=False
        for i in range(len(tract)):
            for j in range(i+1,len(tract)):
                if tract[i,0]>tract[j,0]:
                    save0=tract[i,0]
                    save1=tract[i,1]
                    tract[i,0]=tract[j,0]
                    tract[i,1]=tract[j,1]
                    tract[j,0]=save0
                    tract[j,1]=save1
                    flag=True
    return tract
            