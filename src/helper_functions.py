import numpy as np
import json
from collections import defaultdict
import os, sys
import itertools
import difflib
from glob import glob

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Functions for handling observertions/bed files
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------
def make_callability_from_bed(bedfile, window_size):
    callability = defaultdict(lambda: defaultdict(float))
    with open(bedfile) as data:
        for line in data:

            if not line.startswith('chrom'):

                if len(line.strip().split('\t')) == 3:
                    chrom, start, end = line.strip().split('\t')
                    value = 1
                elif  len(line.strip().split('\t')) > 3:
                    chrom, start, end, value = line.strip().split('\t')[0:4]
                    value = float(value)

                start, end = int(start), int(end)

                firstwindow = start - start % window_size
                lastwindow = end - end % window_size
                
                # not spanning multiple windows (all is added to same window)
                if firstwindow == lastwindow:
                    callability[chrom][firstwindow] += (end-start+1) * value

                # spanning multiple windows
                else:
                    # add to end windows
                    firstwindow_fill = window_size - start % window_size
                    lastwindow_fill = end %window_size
                
                    callability[chrom][firstwindow] += firstwindow_fill * value
                    callability[chrom][lastwindow] += (lastwindow_fill+1) * value

                    # fill in windows in the middle
                    for window_tofil in range(firstwindow + window_size, lastwindow, window_size):
                        callability[chrom][window_tofil] += window_size * value

    return callability


def Load_observations(obs_file, ind, demo_file, window_size = 1000, haploid = True,obs_name="",conditional=False):

    if obs_name=="":
        obs_name=obs_file

    with open(demo_file) as json_file:
        data = json.load(json_file)
    obs_counter = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(list))))
    outgroup_name = []
    state_names = []
    obs= {}
    max_obs=0

    for pop in  (pop for pop in data["pop"] if "ingroup" in pop["type"] ):
        ingroup_name=pop["name"]
    for pop in  (pop for pop in data["pop"] if "outgroup" in pop["type"] ):
        outgroup_name.append(pop["name"])
    for state in  (state for state in data["pop"] if "ancestral" in state["type"] ):
        state_names.append(state["name"])

    
    for pop in  outgroup_name:                     
        with open(f'{obs_file}/{obs_name}.{ind}.{pop}.txt') as dataObs:
            for line in dataObs:
                if not line.startswith('chrom'):
                    chrom, pos, ancestral_base, genotype = line.strip().split()
                    rounded_pos = int(pos)//window_size
                    if haploid:
                        obs_counter[chrom][rounded_pos][pop][0].append(pos)
                        if len(obs_counter[chrom][rounded_pos][pop][0])>max_obs:
                            max_obs=len(obs_counter[chrom][rounded_pos][pop][0])
                    else:
                        for i, base in enumerate(genotype):
                            if base != ancestral_base:
                                obs_counter[chrom][rounded_pos][pop][i].append(pos)
                                if len(obs_counter[chrom][rounded_pos][pop][i])>max_obs:
                                    max_obs=len(obs_counter[chrom][rounded_pos][pop][i])
    i=0
    posMax= {}
    ####
    if conditional:
        times={}
        ancestral={}
        
        for state in state_names:
            ancestries={}
            times = get_split_times(state,data,"outgroup")
            for i in range(len(outgroup_name)):
                ancestries_ingroup = get_ancestries(ingroup_name,data)   
                ancestries_outgtoup = get_ancestries(outgroup_name[i],data)   
                times[i] = max(times[i],get_most_recent_ancestor(ancestries_ingroup,ancestries_outgtoup))
                ancestries[outgroup_name[i]]=times[i]
            ancestries = sorted(ancestries.items(), key=lambda x:x[1], reverse = True)
            ancestral[state]=[]
            for key in ancestries:
                ancestral[state].append(key[0])
        
    ####   
    
    for chrom in obs_counter:
        maxi=-1
        for pos in obs_counter[chrom]:
            if pos>maxi:
                maxi=pos
        posMax[chrom]=maxi  

    if haploid:
        ploidy=1
    else:
        ploidy=2

    if conditional:
        for chrom in obs_counter:
            obs[chrom]=np.zeros((ploidy,posMax[chrom]+1,len(state_names),len(outgroup_name)),dtype=int)
           
            for pos in range(posMax[chrom]+1):         
                for z in range(ploidy):
                    k=0
                    for state in state_names:
                        for i in range(len(outgroup_name)):
                            nb=0
                            for mut in obs_counter[chrom][pos][ancestral[state][0]][z]:
                                b=True
                                for j in range(i+1):
                                    if not mut in obs_counter[chrom][pos][ancestral[state][j]][z]:
                                        b=False
                                if b:
                                    nb+=1
                            obs[chrom][z][pos][k][i]=nb
                        for i in range(len(outgroup_name)-1):
                            obs[chrom][z][pos][k][i]=obs[chrom][z][pos][k][i]-obs[chrom][z][pos][k][i+1]
                        k+=1

    else:
        for chrom in obs_counter:
            obs[chrom]=np.zeros((ploidy,posMax[chrom]+1,len(state_names),len(outgroup_name)),dtype=int)    
            for pos in range(posMax[chrom]+1):         
                for z in range(ploidy):
                    for k in range(len(state_names)):
                        for i in range(len(outgroup_name)):
                            nb=0
                            for mut in obs_counter[chrom][pos][outgroup_name[i]][z]:
                                nb+=1
                            obs[chrom][z][pos][k][i]=nb
                        
    return (obs,max_obs+1)
            

def Load_observations_weights_mutrates(obs_file, weights_file, mutrates_file, window_size = 1000, haploid = False):

    obs_counter = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    haplotypes = defaultdict(int)

    with open(obs_file) as data:
        for line in data:
            if not line.startswith('chrom'):
                chrom, pos, ancestral_base, genotype = line.strip().split()
                rounded_pos = int(pos) - int(pos) % window_size

                if haploid:  
                    for i, base in enumerate(genotype):
                        if base != ancestral_base:
                            obs_counter[chrom][rounded_pos][f'_hap{i+1}'].append(pos)
                            haplotypes[f'_hap{i+1}'] += 1
                else:
                    obs_counter[chrom][rounded_pos][''].append(pos)
                    haplotypes[''] += 1


    chroms, starts, variants, obs = [], [], [], []
    # In the case that there are NO derived variants - use weights to make list of zeros
    if len(obs_counter) == 0:

        if weights_file is None:
            sys.exit(f'{obs_file} is empty! You need to provide a bed file!')

        haplotypes[''] += 1
        callability = make_callability_from_bed(weights_file, window_size)
        for chrom in sorted(callability, key=sortby):
            lastwindow = max(callability[chrom]) + window_size

            for window in range(0, lastwindow, window_size):
                obs_counter[chrom][window][''].append('')
                chroms.append(f'{chrom}')   
                starts.append(window)
                variants.append('')  
                obs.append(0) 

    # Otherwise fill out as normal
    else:
        for haplotype in sorted(haplotypes, key=sortby_haplotype):
            
            for chrom in sorted(obs_counter, key=sortby):
                lastwindow = max(obs_counter[chrom]) + window_size

                for window in range(0, lastwindow, window_size):
                    chroms.append(f'{chrom}{haplotype}')   
                    starts.append(window)
                    variants.append(','.join(obs_counter[chrom][window][haplotype]))  
                    obs.append(len(obs_counter[chrom][window][haplotype])) 
                

    # Read weights file is it exists - else set all weights to 1
    if weights_file is None:
        weights = np.ones(len(obs)) 
    else:  
        callability = make_callability_from_bed(weights_file, window_size)
        weights = []
        for haplotype in sorted(haplotypes, key=sortby_haplotype):
            
            for chrom in sorted(obs_counter, key=sortby):
                lastwindow = max(obs_counter[chrom]) + window_size

                for window in range(0, lastwindow, window_size):
                    weights.append(callability[chrom][window] / float(window_size))


    # Read mutation rate file is it exists - else set all mutation rates to 1
    if mutrates_file is None:
        mutrates = np.ones(len(obs)) 
    else:  
        callability = make_callability_from_bed(mutrates_file, window_size)
        mutrates = []
        for haplotype in sorted(haplotypes, key=sortby_haplotype):
            
            for chrom in sorted(obs_counter, key=sortby):
                lastwindow = max(obs_counter[chrom]) + window_size

                for window in range(0, lastwindow, window_size):
                    mutrates.append(callability[chrom][window] / float(window_size))

    # Make sure there are no places with obs > 0 and 0 in mutation rate or weight
    for index, (observation, w, m) in enumerate(zip(obs, weights, mutrates)):
        if w*m == 0 and observation != 0:
            print(f'warning, you had {observation} observations but no called bases/no mutation rate at index:{index}. weights:{w}, mutrates:{m}')
            obs[index] = 0
            


    return np.array(obs).astype(int), chroms, starts, variants, np.array(mutrates).astype(float), np.array(weights).astype(float)




# ----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Various helper functions
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------
def load_fasta(fasta_file):
    '''
    Read a fasta file with a single chromosome in and return the sequence as a string
    '''
    fasta_sequence = ''
    with open(fasta_file) as data:
        for line in data:
            if not line.startswith('>'):
                fasta_sequence += line.strip().upper()

    return fasta_sequence


def get_consensus(infiles):
    '''
    Find consensus prefix, postfix and value that changes in set of files:

    myfiles = ['chr1.vcf', 'chr2.vcf', 'chr3.vcf']
    prefix, postfix, values = get_consensus(myfiles) 
    
    prefix=chr
    postfix=.vcf
    values = {1,2,3} -> is a set
    '''
    infiles = [str(x) for x in infiles]

    if len(infiles) > 1:

        consensus_strings = defaultdict(int)
        for a, b in itertools.combinations(infiles,2):
            consensus_a = 'START'
            for i,s in enumerate(difflib.ndiff(a, b)):
                if s[0] != ' ':
                    consensus_a += ' '
                else:
                    consensus_a += s[-1]
            consensus_a += 'END'


            new_joined = '|'.join(consensus_a.split()).replace('START','').replace('END','')
            consensus_strings[new_joined] += 1

        for value in consensus_strings:

            if len(value.split('|')) == 2:
                prefix, postfix = value.split('|')
                matches = len([x for x in infiles if prefix in x and postfix in x])

                if matches == len(infiles):
                    values = [x.replace(prefix, '').replace(postfix,'') for x in infiles]
                    return prefix, postfix, set(values)
    else:
        return None, None, None

def sortby(x):
    '''
    This function is used in the sorted() function. It will sort first by numeric values, then strings then other symbols

    Usage:
    mylist = ['1', '12', '2', 3, 'MT', 'Y']
    sortedlist = sorted(mylist, key=sortby)
    returns ['1', '2', 3, '12', 'MT', 'Y']
    '''

    lower_case_letters = 'abcdefghijklmnopqrstuvwxyz'
    if x.isnumeric():
        return int(x)
    elif type(x) == str and len(x) > 0:
        if x[0].lower() in lower_case_letters:
            return 1e6 + lower_case_letters.index(x[0].lower())
        else:
            return 2e6
    else:
        return 3e6
    



def sortby_haplotype(x):
    '''
    This function will sort haplotypes by number
    '''

    if '_hap' in x:
        return int(x.replace('_hap', ''))
    else:
        return x
    



def Make_folder_if_not_exists(path):
    '''
    Check if path exists - otherwise make it;
    '''
    path = os.path.dirname(path)
    if path != '':
        if not os.path.exists(path):
            os.makedirs(path)



def Annotate_with_ref_genome(vcffiles, obsfile):
    obs = defaultdict(list)
    shared_with = defaultdict(str)

    tempobsfile = obsfile + 'temp'

    with open(obsfile) as data, open(tempobsfile,'w') as out:
        for line in data:
            if not line.startswith('chrom'):
                out.write(line)
                chrom, pos, ancestral_base, genotype = line.strip().split()
                derived_variant = genotype.replace(ancestral_base, '')[0]
                ID = f'{chrom}_{pos}'
                obs[ID] = [ancestral_base, derived_variant]

    print('Loading in admixpop snp information')
    for vcffile in handle_infiles(vcffiles):
        command = f'bcftools view -a -R {tempobsfile} {vcffile}'
        print(command)

        for line in os.popen(command):
            if line.startswith('#CHROM'):
                individuals_in_vcffile = line.strip().split()[9:]

            if not line.startswith('#'):

                chrom, pos, _, ref_allele, alt_allele = line.strip().split()[0:5]
                ID =  f'{chrom}_{pos}'
                genotypes = [x.split(':')[0] for x in line.strip().split()[9:]]
                all_bases = [ref_allele] + alt_allele.split(',')

                ancestral_base, derived_base = obs[ID]
                found_in = []

                for original_genotype, individual in zip(genotypes, individuals_in_vcffile):

                    if '.' not in original_genotype:
                        genotype = convert_to_bases(original_genotype, all_bases)   

                        if genotype.count(derived_base) > 0:
                            found_in.append(individual)

                if len(found_in) > 0:
                    shared_with[ID] = '|'.join(found_in)


    # Clean log files generated by vcf and bcf tools
    clean_files('out.log')
    clean_files(tempobsfile)

    return shared_with, individuals_in_vcffile

def handle_individuals_input(argument, group_to_choose):
    if os.path.exists(argument):
        with open(argument) as json_file:
            data = json.load(json_file)
            return data[group_to_choose]
    else:
        return argument.split(',')


# Check which type of input we are dealing with
def handle_infiles(input):
    file_list = glob(input)
    if len(file_list) > 0:
        return file_list
    else:
        if ',' in input:
            return input.split(',')
        else:
            return [input]

# Clean up
def clean_files(filename):
    if os.path.exists(filename):
        os.remove(filename)


# Find variants from admixed population
def flatten_list(variants_list):

    flattened_list = []
    for bin in variants_list:
        if bin != '':
            if ',' in bin:
                for position in bin.split(','):
                    flattened_list.append(position)
            else:
                flattened_list.append(bin)

    return ','.join(flattened_list)


def convert_to_bases(genotype, both_bases):

    return_genotype = 'NN'
    separator = None
    
    if '/' in genotype or '|' in genotype:
        separator = '|' if '|' in genotype else '/'

        base1, base2 = [x for x in genotype.split(separator)]
        if base1.isnumeric() and base2.isnumeric():
            base1, base2 = int(base1), int(base2)

            if both_bases[base1] in ['A','C','G','T'] and both_bases[base2] in ['A','C','G','T']:
                return_genotype = both_bases[base1] + both_bases[base2]

    return return_genotype




# Check which type of input we are dealing with
def combined_files(ancestralfiles, vcffiles):

    # Get ancestral and vcf consensus
    prefix1, postfix1, values1 = get_consensus(vcffiles)
    prefix2, postfix2, values2 = get_consensus(ancestralfiles)

    
    # No ancestral files
    if ancestralfiles == ['']:
        ancestralfiles = [None for _ in vcffiles]
        return ancestralfiles, vcffiles

    # Same length
    elif len(ancestralfiles) == len(vcffiles):
        return ancestralfiles, vcffiles

    # diff lengthts (both longer than 1)       
    elif len(ancestralfiles) > 1 and len(vcffiles) > 1:
        vcffiles = []
        ancestralfiles = []

        for joined in sorted(values1.intersection(values2), key=sortby):
            vcffiles.append(''.join([prefix1, joined, postfix1]))
            ancestralfiles.append(''.join([prefix2, joined, postfix2]))
        return ancestralfiles, vcffiles

    # Many ancestral files only one vcf    
    elif len(ancestralfiles) > 1 and len(vcffiles) == 1:
        ancestralfiles = []
        
        for key in values2:
            if key in vcffiles[0]:
                ancestralfiles.append(''.join([prefix2, key, postfix2]))

        if len(vcffiles) != len(ancestralfiles):
            sys.exit('Could not resolve ancestral files and vcffiles (try comma separated values)')

        return ancestralfiles, vcffiles

    # only one ancestral file and many vcf files
    elif len(ancestralfiles) == 1 and len(vcffiles) > 1:
        vcffiles = []
        
        for key in values1:
            if key in ancestralfiles[0]:
                vcffiles.append(''.join([prefix1, key, postfix1]))

        if len(vcffiles) != len(ancestralfiles):
            sys.exit('Could not resolve ancestral files and vcffiles (try comma separated values)')


        return ancestralfiles, vcffiles
    else:
        sys.exit('Could not resolve ancestral files and vcffiles (try comma separated values)')

#Get for population state, the time it split from each population. Here the split can be an admixture event.
def get_split_times_recombination(state, data, ref):
    times=[]
    ancestries_state = get_ancestries(state,data)
    for outgroup in  (outgroup for outgroup in data["pop"] if ref in outgroup["type"] ):
        ancestries_outgtoup=get_ancestries(outgroup["name"],data)
        times.append(get_most_recent_ancestor(ancestries_state,ancestries_outgtoup))
    
    return times

#Get for population state, the time it split from each population. Here the split can be an admixture even, only for the population which received the admixture.
def get_split_times(state, data, ref):
    times=[]
    for pop in  (pop for pop in data["pop"] if "ingroup" in pop["type"] ):
        ingroup_name=pop["name"]

    ancestries_state = get_ancestries(state,data,state!=ingroup_name)     
    for outgroup in  (outgroup for outgroup in data["pop"] if ref in outgroup["type"] ):
        ancestries_outgtoup=get_ancestries(outgroup["name"],data)
        times.append(get_most_recent_ancestor(ancestries_state,ancestries_outgtoup))
    
    return times

#Get for population state the percentage of its genome which comes from each population.
def get_ancestral_proportion(state, data):
    
    for pop in  (pop for pop in data["pop"] if "ingroup" in pop["type"] ):
        ingroup_name=pop["name"]
        
    ancestral_proportion = {ingroup_name:[0,1]} #At time 0, ingroup as a proportion 1 of ingroup ancestry.
    b=True
    while b:
        b=False
        for elt in data["admixture"]:
            for i in [0,1]:
                pop = elt["ancestral"][i]
                if elt["derived"] in ancestral_proportion and ( not (pop in ancestral_proportion) or (pop in ancestral_proportion and ancestral_proportion[pop][0]<elt["time"] )):
                    if not (pop in ancestral_proportion):
                        ancestral_proportion[pop]=[elt["time"],elt["proportions"][i]*ancestral_proportion[elt["derived"]][1]]
                    else:
                        x=ancestral_proportion[pop][1]
                        y=elt["proportions"][i]
                        ancestral_proportion[pop]=[elt["time"],(y+(1-y)*x)*ancestral_proportion[elt["derived"]][1]]
                    b=True
                    
        for elt in data["split"]:
            for pop in elt["derived"]:
                if pop in ancestral_proportion and ( not (elt["ancestral"] in ancestral_proportion)):
                    ancestral_proportion[elt["ancestral"]]=[elt["time"],ancestral_proportion[pop][1]]
                    b=True

    for elt in ancestral_proportion:
        if elt == state:
            return ancestral_proportion[elt][1]

#Get for population state, each of its ancestral population and the most recent time they interacted with each other (through admixture or split).
def get_ancestries(state,data,admix=True):
    ancestries_state = {state:0}
    b=True
    while b:
        b=False
        for elt in data["admixture"]:
            for i in [0,1]:
                pop = elt["ancestral"][i]
                if admix or elt["proportions"][i]>0.5:
                    if elt["derived"] in ancestries_state and ( not (pop in ancestries_state) or (pop in ancestries_state and ancestries_state[pop]>elt["time"] )):
                        ancestries_state[pop]=elt["time"]
                        b=True          
                    
        for elt in data["split"]:
            for pop in elt["derived"]:
                if pop in ancestries_state and ( not (elt["ancestral"] in ancestries_state) or (elt["ancestral"] in ancestries_state and ancestries_state[elt["ancestral"]]>elt["time"] )):
                    ancestries_state[elt["ancestral"]]=elt["time"]
                    b=True
    
    return ancestries_state

def get_most_recent_ancestor(ancestries_0,ancestries_1):
    res=pow(2,16)
    for key_0 in ancestries_0:
        for key_1 in ancestries_1:
            if key_0==key_1:
                if max(ancestries_0[key_0],ancestries_1[key_1])<res:
                    res=max(ancestries_0[key_0],ancestries_1[key_1])
    return res
        


        