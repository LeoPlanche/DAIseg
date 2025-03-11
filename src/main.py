#!/usr/bin/env python3

import argparse
import numpy as np
import os

from hmm_functions import  DecodeModel, write_HMM_to_file, create_HMM_parameters_from_file, Write_Decoded_output
from bcf_vcf import make_out_group, make_ingroup_obs
#from make_test_data import create_test_data
from make_mutationrate import make_mutation_rate
from helper_functions import Load_observations_weights_mutrates, handle_individuals_input, handle_infiles, combined_files, Load_observations


VERSION = '0.6.9'


def print_script_usage():
    toprint = f'''
Script for identifying introgressed archaic segments (version: {VERSION})

> Turorial:
hmmix make_test_data 
hmmix train  -obs=obs.txt -weights=weights.bed -mutrates=mutrates.bed -param=Initialguesses.json -out=trained.json 
hmmix decode -obs=obs.txt -weights=weights.bed -mutrates=mutrates.bed -param=trained.json


Different modes (you can also see the options for each by writing hmmix make_test_data -h):
> make_test_data        
    -windows            Number of Kb windows to create (defaults to 50,000)
    -nooutfiles         Don't create obs.txt, mutrates.bed, weights.bed, Initialguesses.json (defaults to yes)

> mutation_rate         
    -outgroup           [required] path to variants found in outgroup
    -out                outputfile (defaults to mutationrate.bed)
    -weights            file with callability (defaults to all positions being called)
    -window_size        size of bins (defaults to 1 Mb)

> create_outgroup       
    -ind                [required] ingroup/outgrop list (json file) or comma-separated list e.g. ind1,ind2
    -vcf                [required] path to list of comma-separated vcf/bcf file(s) or wildcard characters e.g. chr*.bcf
    -weights            file with callability (defaults to all positions being called)
    -out                outputfile (defaults to stdout)
    -ancestral          fasta file with ancestral information - comma-separated list or wildcards like vcf argument (default none)
    -refgenome          fasta file with reference genome - comma-separated list or wildcards like vcf argument (default none)

> create_ingroup        
    -ind                [required] ingroup/outgrop list (json file) or comma-separated list e.g. ind1,ind2
    -vcf                [required] path to list of comma-separated vcf/bcf file(s) or wildcard characters e.g. chr*.bcf
    -outgroup           [required] path to variant found in outgroup
    -weights            file with callability (defaults to all positions being called)
    -out                outputfile prefix (default is a file named obs.<ind>.txt where ind is the name of individual in ingroup/outgrop list)
    -ancestral          fasta file with ancestral information - comma-separated list or wildcards like vcf argument (default none)

> train                 
    -obs                [required] file with observation data
    -weights            file with callability (defaults to all positions being called)
    -mutrates           file with mutation rates (default is mutation rate is uniform)
    -param              markov parameters file (default is human/neanderthal like parameters)
    -out                outputfile (default is a file named trained.json)
    -window_size        size of bins (default is 1000 bp)
    -haploid            Change from using diploid data to haploid data (default is diploid)

> decode                
    -obs                [required] directory with observation data
    -ind                individual to decode
    -demo               demographic model
    -weights            file with callability (defaults to all positions being called)
    -mutrates           file with mutation rates (default is mutation rate is uniform)
    -param              markov parameters file (default is human/neanderthal like parameters)
    -out                outputfile prefix <out>.hap1.txt and <out>.hap2.txt if -haploid option is used or <out>.diploid.txt (default is stdout)
    -window_size        size of bins (default is 1000 bp)
    -haploid            Change from using diploid data to haploid data (default is diploid)
    -admixpop ADMIXPOP  Annotate using vcffile with admixing population (default is none)
    -extrainfo          Add variant position for each SNP (default is off)
    '''

    return toprint

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------
def main():

    parser = argparse.ArgumentParser(description=print_script_usage(), formatter_class=argparse.RawTextHelpFormatter)

    subparser = parser.add_subparsers(dest = 'mode')

    # Make outgroup
    outgroup_subparser = subparser.add_parser('create_outgroup', help='Create outgroup information')
    outgroup_subparser.add_argument("-ind",help="[required] ingroup/outgrop list (json file) or comma-separated list e.g. ind1,ind2", type=str, required = True)
    outgroup_subparser.add_argument("-vcf",help="[required] path to list of comma-separated vcf/bcf file(s) or wildcard characters e.g. chr*.bcf", type=str, required = True)
    outgroup_subparser.add_argument("-weights", metavar='',help="file with callability (defaults to all positions being called)")
    outgroup_subparser.add_argument("-out", metavar='',help="outputfile (defaults to stdout)", default = '/dev/stdout')
    outgroup_subparser.add_argument("-ancestral", metavar='',help="fasta file with ancestral information - comma-separated list or wildcards like vcf argument (default none)", default='')
    outgroup_subparser.add_argument("-refgenome", metavar='',help="fasta file with reference genome - comma-separated list or wildcards like vcf argument (default none)", default='')

    # Make mutation rates
    mutation_rate = subparser.add_parser('mutation_rate', help='Estimate mutation rate')
    mutation_rate.add_argument("-outgroup", help="[required] path to variants found in outgroup", type=str, required = True)
    mutation_rate.add_argument("-out", metavar='',help="outputfile (defaults to mutationrate.bed)", default = 'mutationrate.bed')
    mutation_rate.add_argument("-weights", metavar='',help="file with callability (defaults to all positions being called)")
    mutation_rate.add_argument("-window_size", metavar='',help="size of bins (defaults to 1 Mb)", type=int, default = 1000000)

    # Make ingroup observations
    create_obs_subparser = subparser.add_parser('create_ingroup', help='Create ingroup data')
    create_obs_subparser.add_argument("-ind", help="[required] ingroup/outgrop list (json file) or comma-separated list e.g. ind1,ind2", type=str, required = True)
    create_obs_subparser.add_argument("-vcf", help="[required] path to list of comma-separated vcf/bcf file(s) or wildcard characters e.g. chr*.bcf", type=str, required = True)
    create_obs_subparser.add_argument("-outgroup", help="[required] path to variant found in outgroup", type=str, required = True)
    create_obs_subparser.add_argument("-weights", metavar='',help="file with callability (defaults to all positions being called)")
    create_obs_subparser.add_argument("-out", metavar='',help="outputfile prefix (default is a file named obs.<ind>.txt where ind is the name of individual in ingroup/outgrop list)", default = 'obs')
    create_obs_subparser.add_argument("-ancestral", metavar='',help="fasta file with ancestral information - comma-separated list or wildcards like vcf argument (default none)", default='')

    # Decode model
    decode_subparser = subparser.add_parser('decode', help='Decode HMM')
    decode_subparser.add_argument("-obs",help="[required] directory with observation data", type=str, required = True)
    decode_subparser.add_argument("-ind",help="[required] individual", type=str, required = True)
    decode_subparser.add_argument("-demo", metavar='',help="[required] demographic model file", required = True)
    decode_subparser.add_argument("-weights", metavar='',help="file with callability (defaults to all positions being called)")
    decode_subparser.add_argument("-mutrates", metavar='',help="file with mutation rates (default is mutation rate is uniform)")
 
    decode_subparser.add_argument("-out", metavar='',help="outputfile prefix <out>.hap1.txt and <out>.hap2.txt if -haploid option is used or <out>.diploid.txt (default is stdout)", default = '/dev/stdout')
    decode_subparser.add_argument("-window_size", metavar='',help="size of bins (default is 1000 bp)", type=int, default = 1000)
    decode_subparser.add_argument("-haploid",help="Change from using diploid data to haploid data (default is diploid)", action='store_true', default = True)
    decode_subparser.add_argument("-admixpop",help="Annotate using vcffile with admixing population (default is none)")
    decode_subparser.add_argument("-extrainfo",help="Add archaic information on each SNP", action='store_true', default = False)
    decode_subparser.add_argument("-conditional",help="Use conditional probability (default is False)", action='store_true', default = False)

    # all model
    all_subparser = subparser.add_parser('all', help='Run HMM from beginning to end')
    all_subparser.add_argument("-ind",help="[required] ingroup/outgrop list (json file) or comma-separated list e.g. ind1,ind2", type=str, required = True)
    all_subparser.add_argument("-vcfOut",help="[required] path to list of comma-separated vcf/bcf file(s) containing the individuals in the outgroups", type=str, required = True)
    all_subparser.add_argument("-vcfIn",help="[required] path to list of comma-separated vcf/bcf file(s) containing the individuals in the ingroup", type=str, required = True)
    all_subparser.add_argument("-demo",help="demographic model file", type=str, required = True)
    all_subparser.add_argument("-weights", metavar='',help="file with callability (defaults to all positions being called)")
    all_subparser.add_argument("-out", metavar='',help="outputfile (defaults to stdout)", default = '/dev/stdout')
    all_subparser.add_argument("-ancestral", metavar='',help="fasta file with ancestral information - comma-separated list or wildcards like vcf argument (default none)", default='')
    all_subparser.add_argument("-refgenome", metavar='',help="fasta file with reference genome - comma-separated list or wildcards like vcf argument (default none)", default='')
    all_subparser.add_argument("-haploid",help="Change from using diploid data to haploid data (default is diploid)", action='store_true', default = False)
    all_subparser.add_argument("-conditional",help="Use conditional probability (default is False)", action='store_true', default = False)
    
    

    args = parser.parse_args()


    # Decode observations using parameters
    # ------------------------------------------------------------------------------------------------------------
    if args.mode == 'decode':

        #obs, chroms, starts, variants, mutrates, weights  = Load_observations_weights_mutrates(args.obs, args.weights, args.mutrates, args.window_size, args.haploid)
        ingroup_individuals = handle_individuals_input(args.ind,'ingroup')
        os.mkdir(args.out+"/decode")
        hmm_parameters = create_HMM_parameters_from_file(args.demo,conditional=args.conditional)
    
        for individual in ingroup_individuals:
            (obs,max_obs) =  Load_observations(args.out, individual, args.demo, 1000, haploid= args.haploid,obs_name="obs/obs",conditional=args.conditional)
            len_obs = 0
            for key in obs:
                len_obs+=len(obs[key][0])
            

            print('-' * 40)
            print(hmm_parameters)  
            print('> number of windows:', len_obs)
            print('> max observation per window:', max_obs)
            print('> Output prefix is',args.out+"/decode/") 
            print('> Window size is',1000, 'bp') 
            print('> Haploid',args.haploid) 
            print('> Conditional',args.conditional) 
            print('-' * 40)
    
            # Find segments and write output
            segments = DecodeModel(obs, hmm_parameters , max_obs)
            Write_Decoded_output(args.out, segments, args.demo , individual)


    # Create outgroup snps (set of snps to be removed)
    # ------------------------------------------------------------------------------------------------------------
    elif args.mode == 'create_outgroup':

        os.mkdir(args.out)
        # Get list of outgroup individuals
        outgroup_individuals = handle_individuals_input(args.ind, 'outgroup')
    
        # Get a list of vcffiles and ancestral files and intersect them
        vcffiles = handle_infiles(args.vcf)
        ancestralfiles = handle_infiles( args.ancestral)
        refgenomefiles = handle_infiles(args.refgenome)

        ancestralfiles, vcffiles = combined_files(ancestralfiles, vcffiles)
        refgenomefiles, vcffiles = combined_files(refgenomefiles, vcffiles)

        for outgroup in outgroup_individuals:
            print('-' * 40)
            print('> Outgroup individuals:', len(outgroup["ind"]))
            print('> Using vcf and ancestral files')
            for vcffile, ancestralfile, reffile in zip(vcffiles, ancestralfiles, refgenomefiles):
                print('vcffile:',vcffile, 'ancestralfile:',ancestralfile, 'reffile:', reffile)
            print()    
            print('> Callability file:',  args.weights)
            print(f'> Writing output to:', args.out+"/outgroup."+outgroup["name"])
            print('-' * 40)
    
    
            make_out_group(outgroup["ind"], args.weights, vcffiles, args.out+"/outgroup."+outgroup["name"], ancestralfiles, refgenomefiles)


    # Create ingroup observations
    # ------------------------------------------------------------------------------------------------------------
    elif args.mode == 'create_ingroup':
        
        os.makedirs(args.out,exist_ok=True)
        # Get a list of ingroup individuals
        ingroup_individuals = handle_individuals_input(args.ind,'ingroup')

        # Get a list of vcffiles and ancestral files and intersect them
        vcffiles = handle_infiles(args.vcf)
        ancestralfiles = handle_infiles(args.ancestral)
        ancestralfiles, vcffiles  = combined_files(ancestralfiles, vcffiles)

        for filename in os.listdir(args.outgroup):
            print('-' * 40)
            print('> Ingroup individuals:', len(ingroup_individuals))
            print('> Using vcf and ancestral files')
            for vcffile, ancestralfile in zip(vcffiles, ancestralfiles):
                print('vcffile:',vcffile, 'ancestralfile:',ancestralfile)
            print()  
            print('> Using outgroup variants from:', args.outgroup)  
            print('> Callability file:', args.weights)
            print(f'> Writing output to file with prefix: {args.out}.<individual>.txt')
            print('-' * 40)
    
    
            make_ingroup_obs(ingroup_individuals, args.weights, vcffiles, args.out, args.outgroup+"/"+filename, ancestralfiles)


    # Estimate mutation rate
    # ------------------------------------------------------------------------------------------------------------
    elif args.mode == 'mutation_rate':
        os.makedirs(args.out,exist_ok=True)
        for filename in os.listdir(args.outgroup):
            pop=filename.split(sep=".")[1]
            print('-' * 40)
            print(f'> Outgroupfile:', filename)
            print(f'> Outputfile is:', args.out+"/mutationrate."+pop+".bed")
            print(f'> Callability file is:', args.weights)
            print(f'> Window size:', args.window_size)
            print('-' * 40)
    
            make_mutation_rate(args.outgroup+"/"+filename, args.out+"/mutationrate."+pop+".bed", args.weights, args.window_size)
        
    
    # all model
    # ------------------------------------------------------------------------------------------------------------
    elif args.mode == 'all':

        #create outgtoup
        
        os.mkdir(args.out)
        # Get list of outgroup individuals
        outgroup_individuals = handle_individuals_input(args.ind, 'outgroup')
    
        # Get a list of vcffiles and ancestral files and intersect them
        vcffiles = handle_infiles(args.vcfOut)
        ancestralfiles = handle_infiles(args.ancestral)
        refgenomefiles = handle_infiles(args.refgenome)

        ancestralfiles, vcffiles = combined_files(ancestralfiles, vcffiles)
        refgenomefiles, vcffiles = combined_files(refgenomefiles, vcffiles)

        for outgroup in outgroup_individuals:
            print('-' * 40)
            print('> Outgroup individuals:', len(outgroup["ind"]))
            print('> Using vcf and ancestral files')
            for vcffile, ancestralfile, reffile in zip(vcffiles, ancestralfiles, refgenomefiles):
                print('vcffile:',vcffile, 'ancestralfile:',ancestralfile, 'reffile:', reffile)
            print()    
            print('> Callability file:',  args.weights)
            print(f'> Writing output to:', args.out+"/outgroup."+outgroup["name"])
            print('-' * 40)
    
            make_out_group(outgroup["ind"], args.weights, vcffiles, args.out+"/outgroup."+outgroup["name"], ancestralfiles, refgenomefiles)

        #make mutation rate
        
        for filename in (x for x in os.listdir(args.out) if x[0:8]=="outgroup"):
            pop=filename.split(sep=".")[1]
            print('-' * 40)
            print(f'> Outgroupfile:', filename)
            print(f'> Outputfile is:', args.out+"/mutationrate."+pop+".bed")
            print(f'> Callability file is:', args.weights)
            print(f'> Window size:', 100000)
            print('-' * 40)
    
            make_mutation_rate(args.out+"/"+filename, args.out+"/mutationrate."+pop+".bed", args.weights, 100000)
      
        #create ingroup
        os.mkdir(args.out+"/obs")
        ingroup_individuals = handle_individuals_input(args.ind,'ingroup')
        # Get a list of vcffiles and ancestral files and intersect them
        vcffiles = handle_infiles(args.vcfIn)
        ancestralfiles = handle_infiles(args.ancestral)
        ancestralfiles, vcffiles  = combined_files(ancestralfiles, vcffiles)

        for filename in (x for x in os.listdir(args.out) if x[0:8]=="outgroup")   :
            print('-' * 40)
            print('> Ingroup individuals:', len(ingroup_individuals))
            print('> Using vcf and ancestral files')
            for vcffile, ancestralfile in zip(vcffiles, ancestralfiles):
                print('vcffile:',vcffile, 'ancestralfile:',ancestralfile)
            print()  
            print('> Callability file:', args.weights)
            print(f'> Writing output to file with prefix: {args.out}/obs.<individual>.txt')
            print('-' * 40)
    
    
            make_ingroup_obs(ingroup_individuals, args.weights, vcffiles, args.out, args.out+"/"+filename, ancestralfiles,outname="obs/obs")

        #decode
        os.mkdir(args.out+"/decode")
        hmm_parameters = create_HMM_parameters_from_file(args.demo,conditional=args.conditional)

        for individual in ingroup_individuals:
            (obs,max_obs) =  Load_observations(args.out, individual, args.demo, 1000, haploid= args.haploid,obs_name="obs/obs",conditional=args.conditional)
            len_obs = 0
            for key in obs:
                len_obs+=len(obs[key][0])
            
            print('-' * 40)
            print(hmm_parameters)  
            print('> number of windows:', len_obs)
            print('> max observation per window:', max_obs)
            print('> Output prefix is',args.out+"/decode/") 
            print('> Window size is',1000, 'bp') 
            print('> Haploid',args.haploid)
            print('> Conditional',args.conditional)
            print('-' * 40)
    
            # Find segments and write output
            segments = DecodeModel(obs, hmm_parameters , max_obs)
            Write_Decoded_output(args.out, segments, args.demo , individual)

    # Print usage
    # ------------------------------------------------------------------------------------------------------------
    else:
        print(print_script_usage())


if __name__ == "__main__":
    main()