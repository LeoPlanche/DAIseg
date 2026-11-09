# DAIseg
## Overview

HMM-MID is an Hidden Markov Model for Local Ancestry Inference (LAI). It was first designed for the simultaneous detection of multiple archaic introgression in modern humans, but can be used as well to detect any combinaison of modern / archaic introgressions. The user provides a demographic model (defined below) containing a list of outgroups, one ingroup population and the list of the ancestries he is looking for in the ingroup. HMM-MID then uses all the outgroups informations to infer for the ingroup individuals the genomics segments belonging to each of the provided ancestries. An outgroup population may be in the list of ancestries (in that case we would call it a reference population), but an ancestry can be detected even without any reference population, as it is often the case for Archaic detection, our model is agnostic in that regard. The whole model is described in [Planche et al. 2025](https://www.biorxiv.org/content/10.1101/2025.03.17.643330v1) and is based on ideas from [Skov et al. 2018](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1007641), our code is built from [https://github.com/LauritsSkov/Introgression-detection](https://github.com/LauritsSkov/Introgression-detection) where we added the extra layer which allows for the use of demographic models, multiple outgroups and the inference of more than two ancestries.

## Preparing the data

The user needs to provide a demograhic model, a list of samples, the associated bcf files and callability file.

The demography is in a .json file and closely written as in [msprime](https://tskit.dev/msprime/docs/stable/demography.html), it contains:
 - A list of populations, `pop`, each population can take attributes `outgroup`, `ingroup` or `ancestral`. Attributes `outgroup` and `ingroup` are incompatible with each other, but a population may have both attribute `ancestral` and `ingroup` or `outgroup`. To define the non admixed ancestry, it is necessary to add attribute `ancestral` to the population with attribute `ingroup`.
   
   Let us give two quick examples:
   
   To detect Neanderthal, Denisovan and non introgressed segments in Papuans, the `outgroup` populations would be Africans and Eurasians (see [cite]), the `ingroup` would be Papuans and the `ancestral` Neanderthal, Denisovan and Papuans (non introgressed).
   
   To detect Neanderthal and non introgressed segments in modern europeans, the `outgroup` populations could be Africans and Neanderthals (see [Planche et al. 2024](https://link.springer.com/article/10.1134/S1995080224602959)), the `ingroup` modern Europeans and the `ancestral` Neanderthal and modern Europeans.
 - A list of admixtures, `admixture`, containing attributes `time`, the time of admixture, `derived`, the population resulting from admixture, `ancestral`, a list of the two populations which admixed together, and `proportions`, a list of two number, corresponding to the admixture proportion of each of the two ancestral populations.
 - A list of split, `split`, containing attributes `time`, the time of split, `derived`, a list of the two populations resulting from the split, `ancestral`, the population which split into the two derived ones.

Here is a simple example, which corresponds to the detection of Neanderthal in modern Eurasians, using only Africans as outgroup:

```
{
  "pop":[ 
    {"name": "Origin", "type":[]}, 
    {"name": "Sapiens", "type":[]},
    {"name": "Neanderthal", "type":["ancestral"]},
    {"name": "African", "type":["outgroup"]},
    {"name": "Eurasian_before_admixture", "type":[]},
    {"name": "Eurasian", "type":["ingroup","ancestral"]}
  ],
  "admixture":[
    {"time":1850,"derived":"Eurasian","ancestral":["Eurasian_before_admixture", "Neanderthal"], "proportions":[0.97, 0.03]}
  ], 
  "split":[
    {"time":2000,"derived":["African", "Eurasian_before_admixture"],"ancestral":"Sapiens"},
    {"time":30000,"derived":["Sapiens", "Neanderthal"],"ancestral":"Origin"}
  ]
}
```

The list of samples is also given as a .json file. For each population, a list of sample, the names of the population should be the same as in the demographic model.

```
{
 "ingroup":["HG00096","HG00097","HG00099","HG00100","HG00101","HG00102","HG00103"],
"outgroup":[
 {"name":"African", "ind": ["NA18486","NA18488","NA18499","NA18501","NA18502","NA18504","NA18505","NA18507","NA18508"]}
}
```

The variants for the ingroup / outgtoup individuals should be given as a bcf file, the conversion from vcf can be done with bcftools.

The callability file should be given as a .bed file.

## Running the model

Simply run:

```
./main.py all  -demo <demographic>.json -ind <sample list>.json -vcfIn <data for all samples in the ingroup>.bcf -vcfOut <data for all samples in the outgroups>.bcf -weights <mask>.bed -out=<out directory>
```

The bcf files for the ingroup and outgroup can be the same.
The segments can be found for each individuals in:
```
<out directory>/decode/<individual name>.decode.txt.
```

---

If the data is unphased add option `-haploid`.

---

The default decoding algorithm used is Viterbi, add option `-posterior` to use Posterior Decoding.

---

Use `-conditional` to use conditional probabilities as described in Supplementary Materials 1. This option can only be used for specific demographies (see S1), but DAIseg automatically switches back to default if the demography does not allow it, so it can be safely used.



---

Use `-em` to estimate the coalescence times of the demographic model from the data. We describe the algorithm for the specific case of multiple waves of introgression into Papuans in Supplementary Materials 3.

The times are estimated jointly on all the individuals of the ingroup and the segments are then decoded with the estimated values. The `-em` option automatically turns on `-conditional` as well. If the `-conditional` option cannot be used with the demography (see S1), then a warning message is shown and DAIseg automatically switches back to default.

The times of the demographic file are used as the starting point of the estimation. These are the times (given by the user) at which two lineages enter the same population, while the estimated values are the times at which they coalesce, so the estimated values are usually older.  Only the coalescence times used in the emission probabilities are estimated (which corresponds to split times), the admixture times which are used in the transition are not estimated and stay as given in the demographic file.

The algorithm stops after `-em_iter` iterations (default 50), or when the log-likelihood improves by less than `-em_tol` (default 0.1).

The estimated parameters can be found in:
```
<out directory>/em.json
```
and the log-likelihood and the times at each iteration in:
```
<out directory>/em_log.txt
```
