# MHC-peptide-TCR modelling

## 使用说明

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



+ 功能概述: 根据输入的MHC、peptide以及TCR的序列构建复合体结构模型，结构建模方法基于modeller 9.25。

+ 输入项描述：

  -c: 输入的MHC序列所属类型, 包括MHC class I和MHC class II两类可选，输入值必须为<code>MHC-I</code>或 <code>MHC-II</code>；

  -a: MHC分子的alpha chain sequence，可以是包含序列信息的fasta文件，也可以是alpha chain的等位基因名称，输入的基因名称格式如 <code>HLA-A\*01:01, HLA-DPA1\*01:03</code>。程序目前可识别的基因名称见下述两个文件：

  ​	① <code>./program/getTopAlleles/Top_MHCI_alleles.csv</code>文件中存放的为HLA class I每个基因座上的Top10等位基因，由<code>./program/getTopAlleles/getTopMHCIAlleles.py</code>生成。

  ​	② <code>./program/getTopAlleles/Top_MHCII_alleles.csv</code>文件中存放的为HLA class II每个基因座上的Top10等位基因，由<code>./program/getTopAlleles/getTopMHCIIAlleles.py</code>生成。

  ​	**说明：也可以将目前已知的所有等位基因都列举出来供选择**。

  -b: MHC分子的beta chain sequence，可以是包含序列信息的fasta文件，也可以是beta chain的等位基因名称，输入的基因名称格式如 <code>b2m, HLA-DPB1*01:01</code>。程序目前可识别的基因名称见下	述两个文件：

  ​    ① <code>./program/getTopAlleles/b2m.csv</code>文件中存放的为HLA class I中的Beta-2-microglobulin。

  ​	② <code>./program/getTopAlleles/Top_MHCII_alleles.csv</code>文件中存放的为HLA class II每个基因座上的Top10等位基因，由<code>./program/getTopAlleles/getTopMHCIIAlleles.py</code>生成。

  -p: 多肽的序列信息，输入为包含序列信息的fasta文件。

  -A: 可选，TCR的alpha chain sequence，输入为包含序列信息的fasta文件。

  -B: 可选，TCR的beta chain sequence，输入为包含序列信息的fasta文件。

  -u: 用户自定义的输出文件夹。

  -n: 输出多少个模型。

+ 待解决问题：

  TCR目前还没有实现利用基因名称作为输入项。

+ 运行示例 (输出位于<code>model_build</code>文件夹中)：

输入基因名称

```shell
./scripts/Pipeline_modeller_update2021.pl -c MHC-I -a HLA-A*01:01 -b b2m -p ./example/pep_example_sequence.fasta -u MHCI_test -n 1
```



```shell
./scripts/Pipeline_modeller_update2021.pl -c MHC-II -a HLA-DRA*01:01 -b HLA-DRB1*01:01 -p ./example/pep_example_sequence.fasta -u MHCII_test -n 1
```
输入序列文件
```shell
./scripts/Pipeline_modeller_update2021.pl -c MHC-I -a ./example/mhcα_example_sequence.fasta -b ./example/mhcβ_example_sequence.fasta -p ./example/pep_example_sequence.fasta -u MHCI_test_1 -n 1
```
添加TCR序列
```shell
./scripts/Pipeline_modeller_update2021.pl -c MHC-I -a ./example/mhcα_example_sequence.fasta -b ./example/mhcβ_example_sequence.fasta -p ./example/pep_example_sequence.fasta -A ./example/tcrα_example_sequence.fasta -B ./example/tcrβ_example_sequence.fasta -u MHCI_test_2 -n 1
```