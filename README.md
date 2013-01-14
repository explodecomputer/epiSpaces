epiFitness
==========

A genetic algorithm to find 2-locus epistatic patterns that maximise the maintenence of additive genetic variance under selection. This programme was used for analysis in the following article:

*Hemani G, Knott S, Haley C.* **An evolutionary perspective on epistasis and the missing heritability**. *PLoS Genetics (in press)*.

For detailed information on the background and interpretation of `epiFitness` programme please refer to this article.

## Summary

Epistasis may play an important role in the genetic variation of complex traits, as a consequence of being able to mask additive genetic variation from selection. This programme was written to find genotype-phenotype maps that maximise additive genetic variation over many generations, where the causal variants act directly on fitness.

A brief description of what the programme does follows:

We are trying to create two locus bi-allelic epistatic patterns that maximise additive variance and survival time under selection.
- Create `UNMOD` random patterns of epistasis (aka genotype-phenotype maps, or models) with fitness values with a maximum range of URANGE.
- For each model create a population of individuals whose fitness is governed by that model, for each of combinations of frequencies in `UFREQ1` and `UFREQ2`. For example, if `UFREQ1=5` and `UFREQ2=5` then there will be 25 initial populations for each model, and each population will have the mutations at a different starting frequency.
- We calculate the theoretical allele frequency trajectory for each of the `UNFREQS1` x `UNFREQS2` populations, and for each model.
- Additive variance is summed from all the runs of eligable patterns starting from generation `USTARTVA`.
- If at least `UGTHRESH` populations out of all the initial populations survive for at least `UMAXG` generations then that model is eligable to continue.
- Each surviving model is mutated a number of different times (defined in `UNEXTMODS`), by sampling from the values in the `UPERT` array. The original model is also kept so that fitness cannot regress.
- If none survive then choose a new random set of models.


## Installation


Requires GCC. To install (on Mac or Linux) simply clone the repo and run

    make

This will create an execultable called `epiFitness`.


## How to run

`epiFitness` takes only one argument - the filename of a parameter file. An example of a parameter file is available in this repo, `example_parameter.txt`. This algorithm isn't necessarily designed to identify a global solution, moreover it explores the parameter space of the genotype-phenotype map. To this end genetic algorithms can be tuned and perturbed to deliver varying results, and so to access this flexibility the inputs have to be carefully chosen. The input file consists of 19 lines, where each line must be the value for a specific parameter, as listed in order below:


### Parameters

- `UNID` Number of individuals in the population (e.g. **1000**)
- `UMAXG` Maximum number of generations to simulate each population (e.g. **100**)
- `UMAXRUNS` Maximum number of generations for the genetic algorithm to run (e.g. **5000**)
- `UGTHRESH` The number of runs (out of all initial frequencies) that must survive to be allowed to continue (e.g. **20** when there are 25 starting frequencies)
- `USTARTVA` The generation at which to start summing additive genetic variance (e.g. **20**)
- `URANGE` Exact range between maximum and minimum genotype class means after scaling (e.g. **6**)
- `UNFREQ1` Number of different starting frequencies for locus 1 (e.g. **5**)
- `UNFREQ1[0] UNFREQ1[1] ...` Starting frequencies for locus 1 (e.g. **0.1 0.3 0.5 0.7 0.91**)
- `UNFREQ2` Number of different starting frequencies for locus 2 (e.g. **5**)
- `UNFREQ2[0] UNFREQ2[1] ...` Starting frequencies for locus 2 (e.g. **0.1 0.3 0.5 0.7 0.91**)
- `UNMOD` Number of random genotype-phenotype maps to generate (e.g. **40**)
- `UNBESTMODS` Number of models to take forward to the next generation (e.g. **4**)
- `UNEXTMODS[0] UNEXTMODS[1] ... UNEXTMODS[UNBESTMODS-1]` How many random variations of each model carried forward to make, and how many new models to make (e.g. **5 5 5 5 20** if 5 mutations of each of the 4 best models to be made, plus 20 completely new mutations)
- `UNPERT` Number of perturbation values to be used for mutating models (e.g. **10**)
- `UPERT[0] UPERT[1] ...` Perturbation values from which to sample mutations to the models (e.g. **-0.04 -0.02 -0.01 0 0 0 0 0.01 0.02 0.04**)
- `USEED` Seed for random number generator (e.g. **1234**)
- `UFILENAME` Rootname for output files (e.g. **ga_out**)
- `INITIAL` Initial pattern (e.g. **-1 0 1 0 0 0 1 0 -1** would be an example of additive x additive epistasis)
- `OVERRIDE` Overriding pattern for any sampling problems (e.g. **-1 -1 -1 -1 -1 -1 -1 -1 -1**)

Once the input file is created, (e.g. `example_parameters.txt`), the following command runs the programme:

    ./epiFitness example_parameters.txt

### Output files

There are two output files, with the examples given they will be named `ga_out` and `ga_outpats`. `ga_out` is a table that tracks all the models that passed the threshold `UGTHRESH` at each generation, it has the following columns:
- *Generation*
- *Model number*
- *Allele frequency A*
- *Allele frequency B*
- *Average start Vg*
- *Average end Vg*
- *Sum Va over population history*

The `ga_outpats` file details the genotype-phenotype map for every model that appears in `ga_out`, along with its total Va. The patterns that maximise for the conditions of the genetic algorithm are at the end of the files.


## Acknowledgements

This work was produced by [Gibran Hemani][0] part of my PhD thesis at [The Roslin Institute][1] at the University of Edinburgh under the supervision of [Chris Haley][2] and [Sara Knott][3]

## License

EpiFitness is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

EpiFitness is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see [http://www.gnu.org/licenses/][4].

 [0]:http://www.complextraitgenomics.com/the_team/index.php#gibran_hemani
 [1]:http://www.roslin.ac.uk
 [2]:http://www.roslin.ed.ac.uk/chris-haley/
 [3]:http://www.ed.ac.uk/schools-departments/biology/evolutionary-biology/staff-profiles?id=sknott&cw_xml=homepage.php
 [4]:http://www.gnu.org/licenses/

