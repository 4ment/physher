

# Inputs #

## Rules ##
  1. Sequence names in the alignment and tree files must match.
  1. Sequence names cannot contain these characters <font color='red'><code>,[]():;</code></font>
  1. Sequences contain the time of sampling at the end of their names and the date should start with an underscore (e.g. seq1\_2010.5). Sequences are assumed to be isochronous if dates are not present.
  1. No duplicate sequence names.
  1. The tag for calibration points for an internal node inside a tree should be: [&cal\_height={10,15}].
  1. Calibration points are only for isochronous sequences (for now).
  1. NEXUS files:
    * Names can contain spaces but the whole name has to be surrounded by quotes (e.g. "seq 1\_2001.3")
    * Keep the file simple.



## Alignment ##
```
>seq1_2001.3
AAAAAAAA
>seq2_2011.3
AAAAAAAT
>seq3_2013.72
AAAAAAAG
>seq4_2008.2
AAAAAAAC
```

or

```
#NEXUS
begin data;
seq1_2001.3    AATCTCGA
seq2_2011.3    AAACTCGA
seq3_2013.72   AATTTCGA
seq4_2008.2    AATTTCGG
end;
```


## Tree ##
```
(((seq1_2001.3, seq2_2011.3), seq3_2013.72), seq4_2008.2);
```

or

```
#NEXUS
begin trees;
(((seq1_2001.3, seq2_2011.3), seq3_2013.72), seq4_2008.2);
end;
```



# Configuration file #


_physher_ expects a configuration file as argument.

Each line of the configuration file is of the form: key=value

Example:
```
input.sequences = fluA.fa
input.tree      = fluA.tree
output.stem     = fluA
substmodel.type = HKY
clock           = strict
```

# Command line #
Alternatively, _physher_ can be run from the command line:

```
./physher -m HKY -i fluA.fa -t fluA.tree -C strict -o fluA
```
## Inputs/outputs ##

| **Option** | **Key** | **Type** | **Value** | **Default** | **Note** |
|:-----------|:--------|:---------|:----------|:------------|:---------|
| -i | input.sequences | _string_ | path to sequence file | **Mandatory** | sequence file in FASTA or NEXUS format |
| -t | input.tree | _string_ | path to rooted tree file | **Mandatory** |  tree file in newick or NEXUS format |
| -o | output.stem | _string_ | outptut file name | input.sequences |  |

## Substitution model ##

| **Option** | **Key** | **Type** | **Value** | **Default** | **Note** |
|:-----------|:--------|:---------|:----------|:------------|:---------|
| -m | substmodel.type | _string_ | GTR, HKY, K80, JC69, or GY94 | **Mandatory** | 00000, 00001,...,01234 can also be used |
| -r | substmodel.kappa | _real_ | value > 0 | empirical | Only for HKY, K80, and GY94 |
| -r | substmodel.rates | _array of real_ | values > 0 | empirical | Nucleotide only. For GTR it could be:  0.1,0.1,0.2,0.6,0.9 |
|  | substmodel.rates.fixed | _boolean_ | true or false | false | Fixed to values specified in substmodel.rates |
| -f | substmodel.freqs | _array of real_ | 0 < values < 1 | empirical | Nucleotide model only. e.g. 0.1,0.2,0.3,0.4 |
|  | substmodel.freqs.fixed | _boolean_ | true or false| false | Fixed to values specified in substmodel.freqs |
|  | substmodel.codon.omega | _real_ | value > 0  | 1 | GY94 only |
<a href='Hidden comment: 
|| || substmodel.codon.alpha || _real_ || value > 0  || 1 || MG94 only ||
|| || substmodel.codon.beta || _real_ || value > 0  || 1 || MG94 only ||
'></a>
|  | substmodel.heterogeneity | _string_ | no, gamma, inv, or gammainv | no |  |
| -c | substmodel.heterogeneity.gamma.cat | _integer_ | value > 0 | 4 | Number of rate categories |
| -a  | substmodel.heterogeneity.gamma.alpha | _real_ | value > 0  | 0.3 | Parameter of gamma distribution |
|  | substmodel.heterogeneity.gamma.alpha.fixed | _boolean_ | true or false  | false | Fixed to value specified in substmodel.heterogeneity.gamma.alpha |
|  | substmodel.heterogeneity.pinv | _real_ |  0 < value < 1   | 0.5 | proportion of invariant sites |


## Clock ##

| **Option** | **Key** | **Type** | **Value** | **Default** | **Note** |
|:-----------|:--------|:---------|:----------|:------------|:---------|
| -C | clock | _string_ | strict, local, or discrete |  | If not set, no clock is assumed |
| -S | clock.algorithm | _string_ | ga or greedy |  | ga (local and discrete clocks) greedy (local clock only) |
| --forward | clock.forward | _bool_ | true or false | true | Time is forward if a large number = present and smaller number = past |
| --rate | clock.rate | _integer_ | value > 0 |  | It is not fixed, it is just a guess for a better start  |

## Genetic algorithm ##

| **Option** | **Key** | **Type** | **Value** | **Default** | **Note** |
|:-----------|:--------|:---------|:----------|:------------|:---------|
| --ga-pop | ga.popsize | _integer_ | value > 0 | 30 |  |
| --ga-gen  | ga.ngen | _integer_ |  value > 0 | 200 |  |
| --ga-no-improv  | ga.maxnoimprovement | _integer_ | value > 0 | 50 |  |

## Bootstrap ##

| **Option** | **Key** | **Type** | **Value** | **Default** | **Note** |
|:-----------|:--------|:---------|:----------|:------------|:---------|
| -b | bootstrap | _integer_ | value >= 0 | 0 |  |
| --b-threads | bootstrap.threads | _integer_ | value > 0 | 10 |  |

## Miscelleneous ##

| **Option** | **Key** | **Type** | **Value** | **Default** | **Note** |
|:-----------|:--------|:---------|:----------|:------------|:---------|
| -R | random.seed | _integer_ | value >= 0 | random number |  |
| -T | nthreads | _integer_ | value > 0 | 10 | Number of threads used by GA and greedy algorithm|
| --gc | sequences.geneticcode | _integer_ | 0 <= value <= 14 | 0 | Only for codon models. See genetic code [section](Inputfile#Genetic_codes.md) |
| --ic | ic | _string_ | AIC, AICc, or BIC | AICc | For greedy and genetic algorithms |

# Genetic codes #

| **Value** | **Type** |
|:----------|:---------|
| 0 | Universal |
| 1 | Vertebrate Mitochondrial |
| 2 | Yeast |
| 3 | Mold Protozoan Mitochondrial |
| 4 | Mycoplasma |
| 5 | nvertebrate Mitochondrial |
| 6 | Ciliate |
| 7 | Echinoderm Mitochondrial |
| 8 | Euplotid Nuclear |
| 9 | Bacterial |
| 10 | Alternative Yeast |
| 11 | Ascidian Mitochondrial |
| 12 | Flatworm Mitochondrial |
| 13 | Blepharisma Nuclear |
| 14 | No stops |