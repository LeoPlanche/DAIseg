import msprime
import utils
import HMM
import numpy as np
import time
from importlib import reload
from IPython.display import SVG, display


start_time = time.time()

N_DE = 7419
N_AR = 5000
N_AF = 27122
N_Europe = 22000
N_Asia = 24000
N_ME = 4947
N_NE = 7500

Tad = 900 # Admixture time of Denisovans into Papuans
Tsep = 1450 # Split time between Eurasians and Papuans
Tan = 1850 # Admixture time of Neanderthals into Eurasians
Ti = 2000 # Time of Out of Africa
Tn = 24000 # Split time between Denisovan and Neanderthal
Ta = 30000 # Split Archaic and Sapiens

apN = 0.03 #Admixture proportion of Neanderthal into Eurasians
apD = 0.04 #Admixture proportion of Denisovan into Papuans

demography = msprime.Demography()
demography.add_population(name="O", initial_size=N_AR) # Archaic humans
demography.add_population(name="S", initial_size=N_AR) # Homo sapiens
demography.add_population(name="ND", initial_size=N_NE) # Neandertal + Denisovian
demography.add_population(name="N", initial_size=N_NE) # Neandertal
demography.add_population(name="D", initial_size=N_DE) # Denisovian
demography.add_population(name="Af", initial_size=N_AF) # Africans
demography.add_population(name="EAs", initial_size=N_Europe+N_Asia) # Eurasians after Out Of Africa
demography.add_population(name="EAs2", initial_size=N_Europe+N_Asia) # Eurasians after Neanderthal introgression
demography.add_population(name="EAs3", initial_size=N_Europe) # Eurasians after the split with Papuan
demography.add_population(name="P", initial_size=N_ME) # Papuan after the split from Eurasians
demography.add_population(name="P2", initial_size=N_ME) # Papuan after Denisovan introgression

demography.add_admixture(time=Tad, derived="P2", ancestral=["P", "D"], proportions=[1-apD, apD])
demography.add_population_split(time=Tsep, derived=["EAs3", "P"], ancestral="EAs2") 
demography.add_admixture(time=Tan, derived="EAs2", ancestral=["EAs","N"], proportions=[1-apN, apN]) 
demography.add_population_split(time=Ti, derived=["Af", "EAs"], ancestral="S")
demography.add_population_split(time=Tn, derived=["N", "D"], ancestral="ND")
demography.add_population_split(time=Ta, derived=["S", "ND"], ancestral="O")


nbAf = 150
nbEAs3 = 150
nbP2 = 20
seq_len=250000000
rec_rate =1.2e-9
mut_rate=1.25e-8
ploidy=2
L=1000 #Window size for the HMM

# map_file = 'method/recombmap/genetic_map_GRCh37_chr1.txt'
# recomb_map = msprime.RateMap.read_hapmap(map_file)


def createObs(cut,ind,p1,p2,p3) -> list:
    tables = ts.dump_tables()
    nodes = tables.nodes
    seq = np.zeros((int(ts.sequence_length/cut),2),dtype=int) 
    african_id = [p.id for p in ts.populations() if p.metadata['name']=="Af"][0]
    eurasian_id = [p.id for p in ts.populations() if p.metadata['name']=="EAs3"][0] 
    for v in ts.variants():
        i=int(v.site.position/cut)    
        b = c = d =False
        x = 0
        while x < len(v.genotypes):
            if(x==ind):
                if(v.genotypes[x]==1):
                    b = True
                else:
                    x=len(v.genotypes)
            elif(b and nodes[x].population==african_id):
                if(v.genotypes[x]==1):
                    c = True
                    x=p1+p2-1
            elif(b and nodes[x].population==eurasian_id):
                if(v.genotypes[x]==1):
                    d = True 
                    x=len(v.genotypes)
            x=x+1
        if b and not c :
            seq[i][0]=seq[i][0]+1
        if b and not c and not d :
            seq[i][1]=seq[i][1]+1
    return seq

''' 
Initialisation of the start, transition and observation probability matrices
'''
def initS(a,b) -> np.array:
    S = np.zeros(3)
    S[0]=1-a-b
    S[1]=a
    S[2]=b
    return S

def initA(Tn,Td,r,L,a,b) -> np.array:
    A = np.zeros((3,3))
    A[0][1]=Tn*r*L*a
    A[0][2]=Td*r*L*b
    A[1][0]=Tn*r*L*(1-a-b)
    A[1][2]=Td*r*L*b
    A[2][0]=Td*r*L*(1-a-b)
    A[2][1]=Td*r*L*a
    A[0][0]=1-A[0][1]-A[0][2]
    A[1][1]=1-A[1][0]-A[1][2]
    A[2][2]=1-A[2][0]-A[2][1]
    return A

def initB(m,L,Ti,Ta,Tsep,Tn) -> np.array: 
    B = np.empty(shape=(3,80,80))
    meani = m*L*Ti
    meana = m*L*Ta
    meansep = m*L*Tsep
    meann = m*L*Tn
    Psep = np.empty(80)
    Pi_sep = np.empty(80)
    Pa_sep = np.empty(80)
    Pn=np.empty(80)
    Pa_n=np.empty(80)
    Psep[0]=np.exp(-meansep)
    Pi_sep[0]=np.exp(-(meani-meansep))
    Pa_sep[0]=np.exp(-(meana-meansep))
    Pn[0]=np.exp(-(meann))
    Pa_n[0]=np.exp(-(meana-meann))
    sumsep = sumi_sep = suma_sep = sumn = suma_n = 0
    for i in range(1,80):
        Psep[i]=Psep[i-1]*meansep/i
        Pi_sep[i]=Pi_sep[i-1]*(meani-meansep)/i
        Pa_sep[i]=Pa_sep[i-1]*(meana-meansep)/i
        Pn[i]=Pn[i-1]*meann/i
        Pa_n[i]=Pa_n[i-1]*(meana-meann)/i
        sumsep=sumsep+Psep[i]
        sumi_sep=sumi_sep+Pi_sep[i]
        suma_sep=suma_sep+Pa_sep[i]
        sumn=sumn+Pn[i]
        suma_n=suma_n+Pa_n[i]

    Psep[0]=1-sumsep
    Pi_sep[0]=1-sumi_sep
    Pa_sep[0]=1-suma_sep
    Pn[0]=1-sumn
    Pa_n[0]=1-suma_n
  
    for i in range(80):
        for j in range(80):
            if j<=i:
                B[0][i][j]=Psep[j]*Pi_sep[i-j]
                B[1][i][j]=Psep[j]*Pa_sep[i-j]
                B[2][i][j]=Pn[j]*Pa_n[i-j]
            else:
                B[0][i][j]=0
                B[1][i][j]=0
                B[2][i][j]=0
    return B

def printDemography(demography,fileOut):
    with open(fileOut, 'w') as out:
        out.write('{\n')
        out.write('\t"pop": [ \n')
        b=False
        for pop in demography.populations:
            if b:
                out.write(',\n')
            out.write('\t{"name": "'+pop.name+'", "type":['+pop.description+']}')
            b=True

        out.write('\n\t],\n')
        out.write('\t"admixture":[\n')
        b=False
        for event in demography.events:
            if type(event)==msprime.demography.Admixture:
                if b:
                    out.write(',\n')
                b=True
                out.write('\t{"time":'+str(event.time)+',"derived":"'+event.derived+'","ancestral":["'+event.ancestral[0]+'","'+event.ancestral[1]+'"], "proportions":['+str(event.proportions[0])+','+str(event.proportions[1])+']}')
       
        out.write('\n\t],\n')
        out.write('\t"split":[\n')
        b=False
        for event in demography.events:
            if type(event)==msprime.demography.PopulationSplit:
                if b:
                    out.write(',\n')
                b=True
                out.write('\t{"time":'+str(event.time)+',"derived":["'+event.derived[0]+'","'+event.derived[1]+'"],"ancestral":"'+event.ancestral+'"}')
        out.write('\n\t]\n}')

def printIndividuals(demography,fileOut):
    tables = ts.dump_tables()
    nodes = tables.nodes
    with open(fileOut, 'w') as out:
        out.write('{\n')
        out.write('\t"ingroup": [')
        for pop in demography.populations:
            if "ingroup" in pop.description:
                b=False
                for ind in ts.individuals():
                    if nodes[ind.nodes[0]].population == pop.id:
                        if b:
                            out.write(',')
                        b=True
                        out.write('"tsk_'+str(ind.id)+'"')
             
        out.write('],\n')
        out.write('\t"outgroup": [\n')
        c=False
        for pop in demography.populations:
            if "outgroup" in pop.description:
                if c:
                    out.write(',\n')
                c=True
                out.write('\t{"name": "'+pop.name+'", "ind":[')
                b=False
                for ind in ts.individuals():
                    if nodes[ind.nodes[0]].population == pop.id:
                        if b:
                            out.write(',')
                        b=True
                        out.write('"tsk_'+str(ind.id)+'"')
                out.write(']}')

def printMask(seq_len,fileOut):
    with open(fileOut, 'w') as out:
        out.write('1\t0\t'+str(seq_len))
        
M=[]
nbExp=10

for z in range(nbExp):
    ts = msprime.sim_ancestry({"P2": nbP2,"Af":nbAf,"EAs3": nbEAs3 }, 
		                      ploidy=ploidy, 
		                      sequence_length=seq_len,
		                      recombination_rate=rec_rate, 
		                      #recombination_rate = recomb_map,     
		                      demography=demography,record_migrations=True, random_seed=123456789+z)

    ts = msprime.sim_mutations(ts, rate=mut_rate, random_seed=987654321+z)    

    allTractsN = utils.extract_introgressed_tracts(utils.find_introgressed_trees(ts,Tan,'N','P2'))
    allTractsD = utils.extract_introgressed_tracts(utils.find_introgressed_trees(ts,Tad,'D','P2'))
    	
    for w in range(nbP2*ploidy):
		#Get the true tracts from msprime
        tractsD = ((allTractsD[allTractsD['hap_id'] == w])[['start', 'end']] // 1000).astype(int).values.tolist()
        tractsN = ((allTractsN[allTractsN['hap_id'] == w])[['start', 'end']] // 1000).astype(int).values.tolist()
        if len(tractsN)>0 and len(tractsD)>0:
            tractsAf = utils.substract_tracts([[0,seq_len//1000]],np.concatenate((tractsN, tractsD)))
        elif len(tractsN)>0:
            tractsAf = utils.substract_tracts([[0,seq_len//1000]],tractsN)
        elif len(tractsD)>0:
            tractsAf = utils.substract_tracts([[0,seq_len//1000]],tractsD)
        else:
            tractsAf=[[0,seq_len//1000]]
		
        resTrue=[]                                                     
        for i in range(seq_len//1000):
            if utils.inTracts(i,tractsD):
                resTrue.append(2)
            if utils.inTracts(i,tractsN):
                resTrue.append(1)
            if utils.inTracts(i,tractsAf):
                resTrue.append(0)

		#Create observations
        seq=createObs(L,w,nbP2*ploidy,nbAf*ploidy,nbEAs3*ploidy)
		#The states of the HMM, 0: Non archaic, 1: Neanderthal, 2: Denisovan
        states = (0,1,2) 
        cutoff=0.5
        S = initS(apN,apD)
        A = initA(Tan,Tad,rec_rate,L,apN,apD)
        B = initB(mut_rate,L,Ti,Ta,Tsep,Tn)
        resV =  HMM.posterior(seq, S, A,B,cutoff)
        tractsHMM = utils.get_HMM_tracts(resV)
		#Compute the confusion for the current individual
        M.append(utils.confusionMatrix(resV,[tractsAf,tractsN,tractsD]))
       
        
   

