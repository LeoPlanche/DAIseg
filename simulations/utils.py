import numpy as np
import tskit
import pandas as pd


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
            
# Define a function to extract introgressed segments. 
def find_introgressed_trees(
    ts, introgression_time, intro, target,
):    
    intro_population = [p.id for p in ts.populations() if p.metadata['name']==intro][0]
    target_populations = [p.id for p in ts.populations() if p.metadata['name']==target]
    # Extract all information needed for computing introgressed segments from the migration table.
    migration_table        = ts.dump_tables().migrations # The entire migration table.
    migration_nodes        = migration_table.node # Array of ancestral nodes.
    migration_destinations = migration_table.dest # Donor nodes (dest -> source).
    migration_times        = migration_table.time # Recipient nodes (dest -> source).
    migration_left         = migration_table.left # Starting positions [left, right) for migrating segments.
    migration_right        = migration_table.right # Stop positions [left, right) for migrating segments.
    # Extract target haplotypes in the form of list.
    tar_haplotype = []
    for pop in target_populations:
        tar_haplotype.extend(ts.get_samples(pop))
    tar_haplotype = sorted(tar_haplotype)
    # Masks: Introgressing population and introgression time.
    introgression_mask = (migration_destinations == intro_population) # Nodes that belong to the donor population.
    # Donor nodes that correspond to the time of gene flow.
    within_time_mask   = (migration_times[introgression_mask] == introgression_time)
    # Extract the introgressing nodes for the given introgression time and tree coordinates.
    intro_nodes = migration_nodes[introgression_mask][within_time_mask]
    intro_left  = migration_left[introgression_mask][within_time_mask]
    intro_right = migration_right[introgression_mask][within_time_mask]
    # Intialize a dictionary.
    intro_dicc = {
        'hap_id': [],
        't_idx': [], 't_l': [], 't_r': [], 't_span': [],
        'i_idx': [], 'i_l': [], 'i_r': [], 'i_span': [],
        'o_len': [], 'o_perct': [],
    }
    # For every tree with the specified haplotypes...
    for i, tree in enumerate(ts.trees(tracked_samples=tar_haplotype)):
        # For every introgressing node...
        for j, node in enumerate(intro_nodes):
          # If the tree interval is covered, partially or fully, by the migrating node then it is introgressed.
          if ((tree.interval[0] < intro_right[j] and intro_left[j] < tree.interval[1])):
                # Extract the leaves of the introgressing node.
                leaves = set(tree.get_leaves(node))
                # If the introgressed tree has at least one of our target haplotypes...
                if tree.num_tracked_samples(node) > 0:
                    # For every target haplotype.
                    for target in tar_haplotype:
                        # If the haplotype is one of our leaves...
                        if target in leaves:
                            # Extract the coordinates for the tree.
                            t_l, t_r = tree.interval
                            # Calculate the amount and percentage of overlap.
                            o_len = min(t_r, intro_right[j]) - max(t_l, intro_left[j])
                            o_pert = o_len / (t_r - t_l)
                            # Update the dictionary.
                            intro_dicc['hap_id'].append(target)
                            intro_dicc['t_idx'].append(i)
                            intro_dicc['t_l'].append(t_l)
                            intro_dicc['t_r'].append(t_r)
                            intro_dicc['t_span'].append(t_r - t_l)
                            intro_dicc['i_idx'].append(j)
                            intro_dicc['i_l'].append(intro_left[j])
                            intro_dicc['i_r'].append(intro_right[j])
                            intro_dicc['i_span'].append(intro_right[j] - intro_left[j])
                            intro_dicc['o_len'].append(o_len)
                            intro_dicc['o_perct'].append(o_pert)
    # Convert the dictionary to a dataframe.
    intro_df = pd.DataFrame(intro_dicc)
    return intro_df
    
# Define a function to extract contiguous introgressed tracts. Faster than get_migrating_tracts.
def extract_introgressed_tracts(intro_df):
    # Intialize a dictionary.
    tract_dicc = {
        'hap_id': [], 'start': [], 'end': [], 'length': [],
    }
    # For every haplotype with introgressed material.
    for hap_id in intro_df.hap_id.unique():
        # Subset the dataframe.
        hap_df = intro_df[intro_df.hap_id == hap_id]
        # Intialize the tract coordinates.
        start, end = None, None
        # For all trees that overlap introgressed material.
        for _, row in hap_df.iterrows():
            # Extract the coordinates.
            left = row.t_l if row.o_perct == 1 else max(row.t_l, row.i_l)
            right = row.t_r if row.o_perct == 1 else min(row.t_r, row.i_r)
            # If this is the start of a new segment.
            if start is None and end is None:
                # Intialize the start of the segment.
                start = left
                end = right
            # Else-if, this is a continuation of an existing segment.
            elif end == left:
                # Adjust the end coordinate.
                end = right
            # Else the current tract is complete.
            else:
                # Update the dictionary.
                tract_dicc['hap_id'].append(hap_id)
                tract_dicc['start'].append(start)
                tract_dicc['end'].append(end)
                tract_dicc['length'].append(end - start)
                # Initialize the new tract coordinates.
                start, end = left, right
        # Add the last tract if it exists.
        if start is not None and end is not None:
            tract_dicc['hap_id'].append(hap_id)
            tract_dicc['start'].append(start)
            tract_dicc['end'].append(end)
            tract_dicc['length'].append(end - start)
    # Convert the dictionary to a dataframe.
    tract_df = pd.DataFrame(tract_dicc)
    return tract_df