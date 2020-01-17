DRA-DRB-modeller-v2 说明文档：

contact: cy_scu@yeah.net
2020

用途：使用modeller同源建模工具构建MHC分子（以及结合肽）的三维结构
用法：
	1. 准备模板库
		1) "Template_MHCI"和"Template_MHCII"文件夹里的"structure"文件夹中为收集的对应型别的MHC分子的晶体结构，每个晶体结构中均结合有抗原多肽。
		注意："structure"文件夹名称在脚本中被直接使用，不要随意更改

		2) 运行脚本 "ExtractSeq_fromMHCIIComplex-allSeqs.pl structure" 提取出结构文件中的序列信息，存放在"sequence"文件夹中。
		注意："sequence"文件夹名称在脚本中被直接使用，不要随意更改

	2. 准确目标序列文件
		1) 参考example_sequence.fasta文件中的序列格式

	3. 执行基于modeller的同源建模
		1) 运行脚本，例如："./Pipeline_modeller-improve.pl Template_MHCI example_sequence.fasta test 5"
			参数说明： 
				Template_MHCI 为模板（包括结构和序列）所在的文件夹；
				example_sequence.fasta 为需要用于建模的多肽链（输入格式为fasta）
				test 为输出结果存放的自定义文件夹，位于model_build文件夹中。
				5 为构建的模型数目
		
		2) 执行流程：
			读入目标序列文件（fasta格式）
				|
			分别将每一个目标序列文件转换为pir格式（输出对应的.ali文件）
				|
			将目标序列与模板库的模板序列进行比对（NWalign），输出比对打分矩阵（matrix.txt文件）
				|
			使用Hungarian算法寻找目标序列与模板序列之间的最优匹配组合（matrix_assign.txt文件），并计算最优匹配组合的平均得分
				|
			根据平均得分对模板进行排序，选取得分最高的模板作为目标序列的同源建模模板
				|
			按照之前获得的目标序列与模板序列的一一对应关系，分别对对每一对序列进行序列比对（调用脚本"salign.py"）
				|
			将比对好的序列文件合并到一个 .ali 文件中
				|
			执行多条链的同源建模（调用脚本"model-multichain.py"）

		3) 输出文件: 存放在"model_build"文件夹中。
