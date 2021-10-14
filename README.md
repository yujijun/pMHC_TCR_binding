# MHC-peptide-TCR modelling

## Description：
Construct a complex structure model based on the input sequence of MHC, peptide and TCR.

## Requirement:
modeller 9.25

## Usage

```
Usage: Pipeline_modeller_update2021.pl -c <mhc class>
                                -a <mhc alpha chain>
                                -b <mhc beta chain>
                                -p <peptide>
                                [-A] <tcr alpha chain>
                                [-B] <tcr beta chain>
                                -u <userDir>
                                -n <model number>
DESCRIPTION：
      -c  the input MHC class type must be "MHC-I" or "MHC-II"!
      -a  input the sequence file (.fasta) or the allele name (eg:HLA-A*01:01, HLA-DPA1*01:03) of mhc alpha chain.
      -b  input the sequence file (.fasta) or the allele name (eg:b2m, HLA-DPB1*01:01) of mhc beta chain.
      -p   input the sequence file (.fasta) of peptide.
      [-A]  input the sequence file (.fasta) of tcr alpha chain.
      [-B]  input the sequence file (.fasta) of tcr beta chain.
      -u   input the output direction.
      -n   input the output model number.
```


## example
 (output in <code>model_build</code> folder)：

### Take HLA alleles as input:

```shell
./scripts/Pipeline_modeller_update2021.pl -c MHC-I -a HLA-A*01:01 -b b2m -p ./example/pep_example_sequence.fasta -u MHCI_test -n 1
```



```shell
./scripts/Pipeline_modeller_update2021.pl -c MHC-II -a HLA-DRA*01:01 -b HLA-DRB1*01:01 -p ./example/pep_example_sequence.fasta -u MHCII_test -n 1
```

### Take sequence as input：
```shell
./scripts/Pipeline_modeller_update2021.pl -c MHC-I -a ./example/mhcα_example_sequence.fasta -b ./example/mhcβ_example_sequence.fasta -p ./example/pep_example_sequence.fasta -u MHCI_test_1 -n 1
```

### With TCR information
```shell
./scripts/Pipeline_modeller_update2021.pl -c MHC-I -a ./example/mhcα_example_sequence.fasta -b ./example/mhcβ_example_sequence.fasta -p ./example/pep_example_sequence.fasta -A ./example/tcrα_example_sequence.fasta -B ./example/tcrβ_example_sequence.fasta -u MHCI_test_2 -n 1
```

## Evaluation about model
![plot](https://github.com/yujijun/pMHC_TCR_binding/blob/master/Evaluation/example_figure.png)