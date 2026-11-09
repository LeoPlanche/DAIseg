from collections import defaultdict
import numpy as np
from numba import njit
import json
import math

from helper_functions import   Make_folder_if_not_exists,  get_split_times, get_ancestral_proportion, get_most_recent_ancestor, get_most_recent_ancestor_event, get_ancestries, get_split_times_recombination, require_conditional_order

# Shortest coalescence time (generations) an emission class may be given; see DemoStructure.hmm.
MIN_COALESCENCE_TIME = 1.0

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------
# HMM Parameter Class
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------
class HMMParam:
    def __init__(self, state_names, starting_probabilities, transitions, emissions, default_state=0): 
        self.state_names = np.array(state_names)
        self.starting_probabilities = np.array(starting_probabilities)
        self.transitions = np.array(transitions)
        self.emissions = np.array(emissions)
        self.default_state = default_state # index of the ingroup (non admixed) state


    def __str__(self):
        out = f'> state_names = {self.state_names.tolist()}\n'
        out += f'> starting_probabilities = {np.matrix.round(self.starting_probabilities, 3).tolist()}\n'
        out += f'> transitions = {np.matrix.round(self.transitions, 5).tolist()}\n'
        out += f'> emissions = {np.matrix.round(self.emissions*1.25e-8, 3).tolist()}'
        return out

    def __repr__(self):
        return f'{self.__class__.__name__}({self.state_names}, {self.starting_probabilities}, {self.transitions}, {self.emissions})'
        
# Demographic parametrization of the HMM. Each (state, outgroup) emission time is the time at which the two lineages
# first share a population; the event (population, entry time) is one parameter, shared by all pairs that meet there.
# theta = coalescence times of these events (generations), rec_times/props = admixture time and proportion of each state.
class DemoStructure:
    def __init__(self, data, conditional=False):
        self.conditional = conditional
        self.outgroup_names = [pop["name"] for pop in data["pop"] if "outgroup" in pop["type"]]
        self.state_names = [pop["name"] for pop in data["pop"] if "ancestral" in pop["type"]]
        for pop in (pop for pop in data["pop"] if "ingroup" in pop["type"]):
            ingroup_name = pop["name"]
        if ingroup_name not in self.state_names:
            raise ValueError(f'The ingroup population {ingroup_name} must also have type "ancestral"')
        self.default_state = self.state_names.index(ingroup_name)

        ancestries_ingroup = get_ancestries(ingroup_name, data)
        self.param_keys = []
        theta = []
        events = []  # events[state][outgroup] = parameter index
        for state in self.state_names:
            ancestries_state = get_ancestries(state, data, state != ingroup_name)
            events.append([])
            for outgroup in self.outgroup_names:
                ancestries_outgroup = get_ancestries(outgroup, data)
                e = get_most_recent_ancestor_event(ancestries_state, ancestries_outgroup)
                e_in = get_most_recent_ancestor_event(ancestries_ingroup, ancestries_outgroup)
                if e_in[0] > e[0]:
                    e = e_in
                key = f'{e[1]}@{e[0]}'
                if key not in self.param_keys:
                    self.param_keys.append(key)
                    theta.append(e[0])
                events[-1].append(self.param_keys.index(key))
        self.theta = np.array(theta, dtype=float)

        # A[state][k] . theta = emission time of class k (conditional: outgroups sorted by decreasing time, as in Load_observations)
        n_out = len(self.outgroup_names)
        self.A = np.zeros((len(self.state_names), n_out, len(theta)))
        self.order = require_conditional_order(data) if conditional else list(range(n_out))
        for s in range(len(self.state_names)):
            order = self.order
            for k in range(n_out):
                self.A[s, k, events[s][order[k]]] += 1
                if conditional and k < n_out - 1:
                    self.A[s, k, events[s][order[k + 1]]] -= 1

        self.rec_times = np.array(get_split_times_recombination(ingroup_name, data, "ancestral"), dtype=float)
        props = [get_ancestral_proportion(state, data) for state in self.state_names]
        for i in range(len(props)):
            if props[i] == 1:
                props[i] = 2 - sum(props)
        self.props = np.array(props, dtype=float)

    def emission_times(self, theta):
        return np.einsum('skp,p->sk', self.A, theta)

    def hmm(self, theta=None, rec_times=None, props=None, rec_rate=1.2e-8, window_size=1000):
        theta = self.theta if theta is None else theta
        rec_times = self.rec_times if rec_times is None else rec_times
        props = self.props if props is None else props
        n = len(self.state_names)
        # Floor at MIN_COALESCENCE_TIME: with conditional emissions two outgroups of a state can have the same
        # coalescence time, and the class between them would then get rate exactly 0, i.e. probability 0 for any
        # count >= 1 -- a single private variant would make the state impossible in that window. The EM applies the
        # same floor (min_time in m_step_emissions), so the two paths stay consistent.
        emissions = [[max(t, MIN_COALESCENCE_TIME) * window_size for t in row]
                     for row in self.emission_times(theta).tolist()]
        transitions = []
        for i in range(n):
            transitions.append([])
            for j in range(n):
                if i != j and rec_times[i] != 0 and rec_times[j] != 0:
                    transitions[-1].append(window_size*rec_rate*min(rec_times[i], rec_times[j])*props[j])
                elif i != j:
                    transitions[-1].append(window_size*rec_rate*max(rec_times[i], rec_times[j])*props[j])
                else:
                    transitions[-1].append(0)
        for i in range(n):
            transitions[i][i] = 1-sum(transitions[i])
        hmm = HMMParam(state_names=self.state_names,
                       starting_probabilities=list(props),
                       transitions=transitions,
                       emissions=emissions,
                       default_state=self.default_state)
        hmm.param_keys = self.param_keys
        hmm.theta = np.array(theta)
        hmm.rec_times = np.array(rec_times)
        return hmm


# Read HMM parameters from a json file
def create_HMM_parameters_from_file(filename,conditional=False):
    with open(filename) as json_file:
        data = json.load(json_file)
    return DemoStructure(data, conditional).hmm()



# Save HMMParam to a json file
def write_HMM_to_file(hmmparam, outfile):
    data = {key: value.tolist() if isinstance(value, np.ndarray) else value for key, value in vars(hmmparam).items()}
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

def DecodeModel(obs, hmm_parameters, max_obs, posterior_decoding=False):
   
    mut_rate=1.25e-8
    B = initB(hmm_parameters.emissions, mut_rate,max_obs)
    segments={}
    for chrom in obs:
        segments[chrom]={}
        for i in range(len(obs[chrom])):
            if posterior_decoding:
                segments[chrom][i] = posterior(obs[chrom][i], hmm_parameters.starting_probabilities, hmm_parameters.transitions , B, max_obs, cutoff=0.5, default_state=hmm_parameters.default_state) 
            else:
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
                        out.write(f'{chrom}\t{start*window_size}\t{i*window_size}\t{i-start}\t{states_name[curr[i-1]]}\n')
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
 
    return S.astype(int).tolist()


def initialize_matrix(dim1,dim2,value=0):
    F = []
    for i in range(0,dim1):
        F.append([])
        for j in range(0,dim2):
            F[i].append(value)
    return F
     
def posterior(sequence,start,transitions,emissions,max_obs,cutoff=0.5,default_state=0):
    res = initialize_matrix(len(start),len(sequence))
    f=forwardlog(sequence,start,transitions,emissions,max_obs)
    b=backwardlog(sequence,start,transitions,emissions,max_obs)
    for i in range(0,len(start)):
        for j in range(0,len(sequence)):
            res[i][j]=np.exp(f[i][j+1]+b[i][j+1]-f[-1][-1])
            
    # A window is assigned to the most likely state if the probability of not being in the ingroup state is above cutoff
    seq=[]
    for i in range(len(res[0])):
        if 1-res[default_state][i]>cutoff:
            seq.append(int(np.argmax([res[k][i] for k in range(len(start))])))
        else:
            seq.append(default_state)
        
    return seq
            

def forwardlog(sequence,start,transitions,emissions,max_obs):
    M = transitions.shape[0]
    F = initialize_matrix(len(start),len(sequence)+2)
    F[0][0] = 1
    for i in range(0,len(start)):
        F[i][1] = np.log(start[i])+np.log(emissions[i,obs_to_ind(sequence[0][i],M,max_obs)])
    for j in range(2,len(sequence)+1):
        for i in range(0,len(start)):
            p_sum = 0
            sA = []
            for k in range(0,len(start)):
                if  transitions[k,i]*emissions[i,obs_to_ind(sequence[j-1][i],M,max_obs)]==0:
                     p_sum +=0
                else:
                    sA.append(F[k][j-1]+np.log(transitions[k,i])+np.log(emissions[i, obs_to_ind(sequence[j-1][i],M,max_obs)]))
                   
                    
            for  q in range(0,len(sA)):     
                p_sum += np.exp(sA[q]-max(sA))                      
            F[i][j] = np.log(p_sum)+max(sA)
    p_sum = 0
    sA = []
    for k in range(0,len(start)):
        sA.append(F[k][len(sequence)])
    for  q in range(0,len(sA)):     
        p_sum += np.exp(sA[q]-max(sA))  
    F[-1][-1] = np.log(p_sum)+max(sA)
    return F

def backwardlog(sequence,start,transitions,emissions,max_obs):
    M = transitions.shape[0]
    F = initialize_matrix(len(start),len(sequence)+1)
    for i in range(0,len(start)):
        F[i][-1] = 0
    for j in range(len(sequence)-1,0,-1): 
        for i in range(0,len(start)):
            p_sum = 0
            sA=[]
            for k in range(0,len(start)):
                if transitions[i,k]*emissions[k, obs_to_ind(sequence[j][k],M,max_obs)] ==0:
                    p_sum=0
                else:
                    sA.append(F[k][j+1]+np.log(transitions[i,k])+np.log(emissions[k,obs_to_ind(sequence[j][k],M,max_obs)]))
            for  q in range(0,len(sA)):     
                p_sum += np.exp(sA[q]-max(sA))                      
            F[i][j] = np.log(p_sum)+max(sA)      
    p_sum = 0
    sA=[]
    for k in range(0,len(start)):
        if start[k]*emissions[k,obs_to_ind(sequence[0][k],M,max_obs)] == 0:
            p_sum +=0
        else:
            sA.append(F[k][1]+np.log(start[k]*emissions[k,obs_to_ind(sequence[0][k],M,max_obs)]))
    for  q in range(0,len(sA)):     
         p_sum += np.exp(sA[q]-max(sA)) 
    F[0][0] = np.log(p_sum)+max(sA)
    return F
