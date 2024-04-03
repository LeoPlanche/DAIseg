# HHM-MID
## Overview

HMM-MID is an Hidden Markov Model for Local Ancestry Introgression (LAI). It was first designed for the detection of multiple archaic introgression in modern humans, but can be used to detect any combinaison of modern / archaic introgressions as well. The user provides a demographic model, a list of outgrouŝ, a list of reference populations, one ingroup population and the associated vcf files. HMM-MID then uses all the outgroups informations to infer for the ingroup individuals the genomics segments belonging to each of the provided ancestries. An outgroup population may be in the list of ancestries (in that case we would call it a reference population), but an ancestry can be detected even without any reference population, as it is often the case for Archaic detection, our model is agnostic in that regard. The whole model is described [cite] and is based on ideas from [Skov et al. 2018](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1007641), our code is built from [https://github.com/LauritsSkov/Introgression-detection](https://github.com/LauritsSkov/Introgression-detection) where we added the extra layer which allows for the use of demographic models, multiple outgroups and the inference of more than two ancestries.

## Preparing the data

The user needs to provide a demograhic model, a list of samples, the associated bcf files and callability file.

The demography is in a .json file and closely written as in [msprime](https://tskit.dev/msprime/docs/stable/demography.html), it contains:
 - A list of populations, 'pop', each population can take attributes 'outgroup', 'ingroup' or 'ancestral'. Attributes 'outgroup' and 'ingroup' are imcompatible with each other, but a population may have both attribute 'ancestral' and 'ingroup' or 'outgroup'. Let us give two quick examples, to detect Denisovan in Papuans, the outgroup populations would be Africans and Eurasians, the ingroup would be Papuans and the ancestral Neanderthal, Denisovan, non introgressed. To detect Neanderthal in modern europeans, the outgroup populations could be Africans and Neanderthals (see [citeAnna]), the ingroup modern Europeans and the ancestral Neanderthal and non introgressed.
 - A list of admixtures, 'admixture', containing attributes 'time', the time of admixture, 'derived', the population resulting from admixture, 'ancestral', a list of the two populations which admixed together, and 'proportions', a list of two number, corresponding to the admixture proportion of each of the two ancestral populations.
 - A list of split, 'split', containing attributes 'time', the time of split, 'derived', a list of the two populations resulting from the split, 'ancestral', the population which split into the two derived ones.

Here is a simple example, which corresponds to the detection of Neanderthal in modern Eurasians, using only Africans as outgroup:

```
{
  "pop":[ 
    {"name": "Origin", "type":[]}, 
    {"name": "Sapiens", "type":[]},
    {"name": "Neanderthal", "type":["ancestral"]},
    {"name": "African", "type":["outgroup"]},
    {"name": "Eurasian_before_admixture", "type":[]}
    {"name": "Eurasian", "type":["ingroup","ancestral"]}
  ],
  "admixture":[
    {"time":1850,"derived":"Eurasian","ancestral":["Eurasian_before_admixture", "Neanderthal"], "proportions":[0.97, 0.03]}
  ], 
  "split":[
    {"time":2000,"derived":["African", "Eurasian_before_admixture"],"ancestral":"Sapiens"},
    {"time":30000,"derived":["Sapiens", "Neanderthal"],"ancestral":"O"}
  ]
}
```

The list of samples is also given as a .json file. For each population, a list of sample, the names of the population should be the same as in the demographic model.

```
{
 "ingroup":[
  {"name": "Eurasian", "ind": ["HG00096","HG00097","HG00099","HG00100","HG00101","HG00102","HG00103"]}
 ],
"outgroup":[
 {"name":"African", "ind": ["NA18486","NA18488","NA18499","NA18501","NA18502","NA18504","NA18505","NA18507","NA18508"]}
}
```

The variants for the ingroup / outgtoup individuals should be given as a bcf file, the conversion from vcf can be done with bcftools.

The callability file should be given as a .bed file.

## Running the model

Simply run:

```
./main.py whole  -demo=<demographic>.json -ind=<sample list>.json -vcfIn=<data for all samples in the ingroup>.bcf -vcfOut=<data for all samples in the outgroups>.bcf -weights=<mask>.bed -out=<out directory>
```
The bcf files for the ingroup and outgroup can be the same.
The segments can be found for each individuals in:
```
<out directory>/decode/<individual name>.decode.txt.
```
