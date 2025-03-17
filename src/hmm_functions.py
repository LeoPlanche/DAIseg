from collections import defaultdict
import numpy as np
from numba import njit
import json
import math

from helper_functions import  Annotate_with_ref_genome, Make_folder_if_not_exists, flatten_list, get_split_times, get_ancestral_proportion, get_most_recent_ancestor, get_ancestries, get_split_times_recombination

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------
# HMM Parameter Class
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------
class HMMParam:
    def __init__(self, state_names, starting_probabilities, transitions, emissions): 
        self.state_names = np.array(state_names)
        self.starting_probabilities = np.array(starting_probabilities)
        self.transitions = np.array(transitions)
        self.emissions = np.array(emissions)


    def __str__(self):
        out = f'> state_names = {self.state_names.tolist()}\n'
        out += f'> starting_probabilities = {np.matrix.round(self.starting_probabilities, 3).tolist()}\n'
        out += f'> transitions = {np.matrix.round(self.transitions, 5).tolist()}\n'
        out += f'> emissions = {np.matrix.round(self.emissions*1.25e-8, 3).tolist()}'
        return out

    def __repr__(self):
        return f'{self.__class__.__name__}({self.state_names}, {self.starting_probabilities}, {self.transitions}, {self.emissions})'
        
# Read HMM parameters from a json file
def create_HMM_parameters_from_file(filename,conditional=False):
    rec_rate =1.2e-9

    with open(filename) as json_file:
        data = json.load(json_file)

    state_names = []
    outgroup_name = []
    emissions = []
    starting_probabilities = []
    transitions = []
    for pop in  (pop for pop in data["pop"] if "ingroup" in pop["type"] ):
        ingroup_name=pop["name"]

    for pop in  (pop for pop in data["pop"] if "outgroup" in pop["type"] ):
        outgroup_name.append(pop["name"])
        

    for state in  (state for state in data["pop"] if "ancestral" in state["type"] ):
        state_names.append(state["name"])

    if conditional:
        for state in state_names:
            starting_probabilities.append(get_ancestral_proportion(state,data))        
            emissions.append([])
            times = get_split_times(state,data,"outgroup")
            for i in range(len(outgroup_name)):
                ancestries_ingroup = get_ancestries(ingroup_name,data)   
                ancestries_outgtoup = get_ancestries(outgroup_name[i],data)   
                times[i] = max(times[i],get_most_recent_ancestor(ancestries_ingroup,ancestries_outgtoup))
            
            times.sort(reverse=True)
            for i in range(len(times)-1):
                emissions[-1].append((times[i]-times[i+1])*1000)
            emissions[-1].append(times[-1]*1000)
            
        for i in range(len(starting_probabilities)):
            if starting_probabilities[i]==1:
                starting_probabilities[i]=2-sum(starting_probabilities)
    else:
        for state in state_names:
            starting_probabilities.append(get_ancestral_proportion(state,data))        
            emissions.append([])
            times = get_split_times(state,data,"outgroup")
            for i in range(len(outgroup_name)):
                ancestries_ingroup = get_ancestries(ingroup_name,data)   
                ancestries_outgtoup = get_ancestries(outgroup_name[i],data)   
                times[i] = max(times[i],get_most_recent_ancestor(ancestries_ingroup,ancestries_outgtoup))
                emissions[-1].append(times[i]*1000)
            
        for i in range(len(starting_probabilities)):
            if starting_probabilities[i]==1:
                starting_probabilities[i]=2-sum(starting_probabilities)


    recombination_time=get_split_times_recombination(ingroup_name,data,"ancestral")
    for i in range(len(state_names)):
        transitions.append([])
        for j in range(len(state_names)):
            if i!=j and recombination_time[i]!=0 and recombination_time[j]!=0:
                transitions[-1].append(1000*rec_rate*min(recombination_time[i],recombination_time[j])*starting_probabilities[j])
            elif i!=j:
                transitions[-1].append(1000*rec_rate*max(recombination_time[i],recombination_time[j])*starting_probabilities[j])
            else:
                transitions[-1].append(0)
                
    for i in range(len(state_names)):
        transitions[i][i]=1-sum(transitions[i])
        
    return HMMParam(state_names=state_names,
                    starting_probabilities=starting_probabilities,
                    transitions=transitions,
                    emissions=emissions)



# Save HMMParam to a json file
def write_HMM_to_file(hmmparam, outfile):
    data = {key: value.tolist() for key, value in vars(hmmparam).items()}
    json_string = json.dumps(data, indent = 2) 
    with open(outfile, 'w') as out:
        out.write(json_string)


def logoutput(hmm_parameters, loglikelihood, iteration):

    n_states = len(hmm_parameters.emissions)

    # Make header
    if iteration == 0:    
        print_emissions = '\t'.join(['emis{0}'.format(x + 1) for x in range(n_states)])
        print_starting_probabilities = '\t'.join(['start{0}'.format(x + 1) for x in range(n_states)])
        print_transitions = '\t'.join(['trans{0}_{0}'.format(x + 1) for x in range(n_states)])
        print('iteration', 'loglikelihood', print_starting_probabilities, print_emissions, print_transitions, sep = '\t')

    # Print parameters
    print_emissions = '\t'.join([str(x) for x in np.matrix.round(hmm_parameters.emissions, 4)])
    print_starting_probabilities = '\t'.join([str(x) for x in np.matrix.round(hmm_parameters.starting_probabilities, 3)])
    print_transitions = '\t'.join([str(x) for x in np.matrix.round(hmm_parameters.transitions, 4).diagonal()])
    print(iteration, round(loglikelihood, 4), print_starting_probabilities, print_emissions, print_transitions, sep = '\t')






# ----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Decode
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------

def DecodeModel(obs, hmm_parameters, max_obs):
   
    mut_rate=1.25e-8
    B = initB(hmm_parameters.emissions, mut_rate,max_obs)
    segments={}
    for chrom in obs:
        segments[chrom]={}
        for i in range(len(obs[chrom])):
            segments[chrom][i] = viterbi(obs[chrom][i], hmm_parameters.starting_probabilities, hmm_parameters.transitions , B, max_obs) 
    return segments



def Write_Decoded_output(outputprefix, segments, filename , ind="" , window_size=1000):
    with open(filename) as json_file:
        data = json.load(json_file)

    states_name =[]
    for pop in  (pop for pop in data["pop"] if "ancestral" in pop["type"] ):
        states_name.append(pop["name"])

    if ind!="":
        filename=f'{outputprefix}/decode/{ind}.decode.txt'
    else:
        filename=f'{outputprefix}.decode.txt'
        
    with open(filename, 'w') as out:
        out.write('chrom\tstart\tend\tlength\tstate\n')
    
        for chrom in segments:
            for ploidy in segments[chrom]:
                curr=segments[chrom][ploidy]
                start=0
                for i in range(1,len(curr)):
                    if curr[i]!=curr[i-1]:
                        out.write(f'{chrom}\t{start*window_size}\t{(i-1)*window_size}\t{i-start}\t{states_name[curr[i-1]]}\n')
                        start=i
                out.write(f'{chrom}\t{start*window_size}\t{len(curr)*window_size}\t{len(curr)-start}\t{states_name[curr[-1]]}\n')

    out.close()


def initB(emissions, mut_rate, max_obs) -> np.array: 

    nb_states = len(emissions)
    nb_outgroup = len(emissions[0])
    B = np.empty(shape=(nb_states,pow(max_obs,nb_outgroup)))
    B.fill(1)

    proba = np.empty(shape=(nb_states,nb_outgroup,max_obs))
    for i in range(nb_states):
        for j in range(nb_outgroup):
            proba[i][j][0] = np.exp(-emissions[i][j]*mut_rate)
            
    for i in range(nb_states):
        for j in range(nb_outgroup):      
            for k in range(1,max_obs):
                proba[i][j][k]=proba[i][j][k-1]*(emissions[i][j]*mut_rate)/k
                
    for i in range(nb_states):
        for j in range(nb_outgroup):
            proba[i][j][0] = 1-sum(proba[i][j])+proba[i][j][0]

    for i in range(nb_states):
        curr =  np.zeros(shape=(nb_outgroup),dtype=int)
        for j in range(pow(max_obs,nb_outgroup)):
            for k in range(nb_outgroup):
                B[i][j]=B[i][j]*proba[i][k][curr[k]]
            curr=incr(curr,nb_outgroup,max_obs)          
    return B

def incr(array,size,max):
    b=True
    i=0
    while(i<size):
        if array[i]<max-1:
            array[i]+=1
            return array
        else:
            array[i]=0
            i+=1

def obs_to_ind(V,nb_outgroup,max_obs):
    curr=1
    res=0
    for elt in V:
        res+=elt*curr
        curr=curr*max_obs
    return res


def viterbi(V, initial_distribution, a, b, max_obs):
    T = len(V)
    M = a.shape[0]

    
    omega = np.zeros((T, M))
    for i in range(M):
        omega[0, i] = np.log(initial_distribution[i] * b[i][ obs_to_ind(V[0][i],M,max_obs)])
 
    prev = np.zeros((T - 1, M))
    
 
    for t in range(1, T):
        for j in range(M):
            # Same as Forward Probability
            probability=[]
            for i in range(M):
                probability.append(omega[t - 1][i] + np.log(a[i, j]) + np.log(b[j][ obs_to_ind(V[t][j],M,max_obs)]))

 
            # This is our most probable state given previous state at time t (1)
            prev[t - 1, j] = np.argmax(probability)
 
            # This is the probability of the most probable state (2)
            omega[t, j] = np.max(probability)
 
    # Path Array
    S = np.zeros(T)
 
    # Find the most probable last hidden state
    last_state = np.argmax(omega[T - 1, :])
 
    S[0] = last_state
 
    backtrack_index = 1
    for i in range(T - 2, -1, -1):
        S[backtrack_index] = prev[i, int(last_state)]
        last_state = prev[i, int(last_state)]
        backtrack_index += 1
 
    # Flip the path array since we were backtracking
    S = np.flip(S, axis=0)
 
    # Convert numeric values to actual hidden states
 
    result = []
    for s in S:
        if s == 0:
            result.append(0)
        elif s == 1:
            result.append(1)
        elif s == 2:
            result.append(2)
 
    return result

        
    