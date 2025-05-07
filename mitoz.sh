# 🌼🌼🌼 mitogenome组装简单版 🌼🌼🌼
source activate getorganelle
cd /data01/xush/1.mitogenome/0.data
for i in *_1.fq.gz
do
i=${i%_1.fq.gz};
extractfq -fq1 ${i}_1.fq.gz \
          -fq2 ${i}_2.fq.gz \
	  -outfq1 ${i}_1_extractfq.fq.gz \ 
          -outfq2 ${i}_2_extractfq.fq.gz \
	  -size_required 2 -gz
done
---------------------------------------------------
安装mitoz的环境：先下载压缩包，然后激活环境
mkdir -p ./mitoz3.6
tar -xzf ./mitoz3.6.tar.gz -C ./mitoz3.6
下载地址：https://www.dropbox.com/scl/fo/4md8irodd9flywxhp85wf/AEXYFykNDPgVuJ6nSS6CCjs?rlkey=gpouy6vue1rf9dgva4jqqcik0&e=1&dl=0 #直接下载
source /data01/xush/1.mitogenome/mitoz3.6/bin/activate
conda-unpack
---------------------------------------------------
source /data01/xush/1.mitogenome/mitoz3.6/bin/activate
cd /data01/xush/1.mitogenome/0.data
for i in ../0.data/*_1_extractfq.fq.gz; do
    filename=$(basename "$i")              # 取出文件名部分
    sample=${filename%_1_extractfq.fq.gz}  # 去掉后缀变成样本名
    mitoz all \
    --outprefix ${sample} \
    --thread_number 20 \
    --clade Chordata \
    --genetic_code 2 \
    --species_name "Sinocyclocheilus" \
    --fq1 ../0.data/${sample}_1_extractfq.fq.gz \
    --fq2 ../0.data/${sample}_2_extractfq.fq.gz \
    --fastq_read_length 150 \
    --assembler megahit \
    --requiring_taxa Cyprinidae \
    --workdir /data01/xush/1.mitogenome/2.results/${sample}_workdir
done
#或者不通过上面的直接加入参数--data_size_for_mt_assembly 3 
#✅如果是脊椎动物（如鱼类、哺乳类、鸟类等）：
--genetic_code 2
代表 Vertebrate Mitochondrial Code（脊椎动物线粒体遗传密码）
✅如果是无脊椎动物：
昆虫：--genetic_code 5 （Invertebrate Mitochondrial）
其他具体分类可能不同，可以查：NCBI genetic code tables
---------------------------------------------


---------------------------------------------
# 🌼🌼🌼mitogenome组装复杂版🌼🌼🌼
# 质控，fastqc *.fq.gz，cleandata的质量较好
# 用bwa进行Mapping 基因组选择近缘物种的基因组
cd /data01/xush/1.mitogenome/1.mapping/ref
bwa index Sang-hap1.chr.fasta

cd /data01/xush/1.mitogenome/0.data
input=/data01/xush/1.mitogenome/0.data
outdir=/data01/xush/1.mitogenome/1.mapping
ref=/data01/xush/1.mitogenome/1.mapping/ref/Sang-hap1.chr.fasta

for i in *_1.fq.gz; do
    sample=${i%_1.fq.gz}
    fq1=${sample}_1.fq.gz
    fq2=${sample}_2.fq.gz
    bwa mem -t 10 -M -R "@RG\tID:$sample\tPL:illumina\tLB:$sample\tSM:$sample" \
        $ref $input/$fq1 $input/$fq2 | \
    samtools view -Sb -o $outdir/${sample}.bam
done
----------------------------------------------------------------------
## 线粒体组装注释流程 ##
----------------------------------------------------------------------
# Sang-hap1.chr.fasta中没有线粒体序列
# mapping得到bam文件，把其中unmap的部分提取出来：（4+8=12，两端都没比对上的）
# 然后按名称排序，然后用bedtools将bam转为fastq
	#!/bin/bash
	cd /work/users/wangyw/02_S.micro/output/bam
	for i in `ls *bam`;
	do
	samtools view -bf 12 $i > ${i%%.*}-unmap.bam
	samtools sort -n ${i%%.*}-unmap.bam -o ${i%%.*}-unmap-namesort.bam
	bedtools bamtofastq -i ${i%%.*}-unmap2-namesort.bam -fq ../mt_fq/${i%%.*}-unmap_r1.fq -fq2 ../mt_fq/${i%%.*}-unmap_r2.fq
	rm ${i%%.*}-unmap.bam
	rm ${i%%.*}-unmap-namesort.bam
	done
----------------------------------------------------------------------
#### 把unmap的fq序列再map到线粒体基因组上(从网上下载小眼金线鲃线粒体基因组 GenBank: MN145877.1)
	#!/bin/bash
	mydir=/work/users/wangyw/02_S.micro
	genome=$mydir/input/mt/S.micro.fa
	index=$mydir/input/mt/S.micro.fa

	cd $mydir/output/mt_fq

	for i in `ls *unmap_r1.fq`
	do
	fq1=$i
	fq2=${i%-*}-unmap_r2.fq
	sample=${i%-*}

	bwa mem -M -R "@RG\tID:$sample\tPL:illumina\tLB:$sample\tSM:$sample" $index $fq1 $fq2 | samtools view -Sbh -o $mydir/output/mt_bam/$sample.bam
	samtools flagstat $mydir/output/mt_bam/$sample.bam > $mydir/output/mt_bam/${sample}.flagstat
	done
#统计
for i  in `ls *.flagstat`;do reads=`awk 'NR==5' $i` ;echo $i $reads >>mt-raw-stat;done
#删除unmap.fq
----------------------------------------------------------------------
#### 提取线粒体序列：
	#先过滤:
	cd /work/users/wangyw/02_S.micro/output/mt_bam
	for i in `ls *.bam`
	do 
	samtools view -b -f 2 -F 256 -q 30 $i -o ${i%%.*}-f2F256q30.bam
	samtools flagstat ${i%%.*}-f2F256q30.bam > ../mt_bamflagstat/${i%%.*}-f2F256q30.flagstat
	done
	#统计
	for i  in `ls *-f2F256q30.flagstat`;do reads=`awk 'NR==5' $i` ;echo $i $reads >>mt-filt-stat;done
	# 删除*.bam

	#提取线粒体序列，方法同上。
	cd /work/users/wangyw/02_S.micro/output/mt_bam
	for i in `ls *f2F256q30.bam`;
	do
	samtools sort -n $i -o ${i%%.*}-namesort.bam
	bedtools bamtofastq -i ${i%%.*}-namesort.bam -fq ../mt_fq/${i%%.*}-mt_r1.fq -fq2 ../mt_fq/${i%%.*}-mt_r2.fq
	rm ${i%%.*}-namesort.bam
	done
----------------------------------------------------------------------
#### 组装线粒体基因组：
#使用MitoZ来组装（https://github.com/linzhi2013/MitoZ）
	#组装
	(mitozEnv) [wangyw@Buffalo mt_mitoz]$ cat mt_mitoz_all.sh
	mydir=/work/users/wangyw/02_S.micro
	fq1=$1
	fq2=$2
	outfile=$3

	python3 /work/users/wangyw/biosoft/release_MitoZ_v2.3/MitoZ.py all2 \
	--genetic_code 2 \
	--clade Chordata \
	--insert_size 350 \
	--thread_number 8 \
	--fastq1 $fq1 \
	--fastq2 $fq2 \
	--outprefix $outfile \
	--fastq_read_length 150 \
	--run_mode 2 \
	1>${outfile}.log 2>${outfile}.err

	if [[ ! -d "$tmp" ]];then mv tmp ${outfile}.tmp ;fi

	#批量组装脚本
	(mitozEnv) [wangyw@Buffalo mt_mitoz]$ cat all.sh
	for i in `ls /work/users/wangyw/02_S.micro/output/mt_fq/*-f2F256q30-mt_r1.fq`
	do
	filen1=${i%-*};
	filen2=${filen1##*/}
	filen3=${filen2%-*}
	sh mt_mitoz_all.sh ${i%_*}_r1.fq ${i%_*}_r2.fq $filen3
	done

#先进入环境，然后开始组装
source activate mitozEnv
nohup sh all.sh > all.sh.nohup 2>&1 &


##SCP-10没有结果，因为数据较少? 
#换用小的kmer来组装
source activate mitozEnv
(mitozEnv) [wangyw@Buffalo SCP-10.assembly]$ pwd
/work/users/wangyw/02_S.micro/output/tmp/SCP-10.assembly
(mitozEnv) [wangyw@Buffalo SCP-10.assembly]$ /work/users/wangyw/biosoft/release_MitoZ_v2.3/bin/assemble/mitoAssemble all -K 51 -o work51 -s work71.soaptrans.lib -p 8
#提取出mito
(mitozEnv) [wangyw@Buffalo SCP-10.assembly]$ cat findmitoscaf.sh
mydir=/work/users/wangyw/02_S.micro
fq1=$mydir/output/mt_fq/SCP-10-f2F256q30-mt_r1.fq
fq2=$mydir/output/mt_fq/SCP-10-f2F256q30-mt_r2.fq
python3 /work/users/wangyw/biosoft/release_MitoZ_v2.3/MitoZ.py findmitoscaf --genetic_code 2 --clade Chordata \
--outprefix test \
--thread_number 8 \
--fastq1 $fq1 \
--fastq2 $fq2 \
--fastq_read_length 150 \
--fastafile work51.scafSeq
#注释
python3 /work/users/wangyw/biosoft/release_MitoZ_v2.3/MitoZ.py annotate --genetic_code 2 --clade Chordata \
--outprefix SCP-10-test --thread_number 8 \
--fastafile test.mitogenome.fa
#最后的结果文件为SCP-10-test.fasta
----------------------------------------------------------------------

#使用MitoS来注释
#网址：http://mitos2.bioinf.uni-leipzig.de/index.py

#一些需要进行反向互补
source activate mitozEnv
for i in `ls *fasta`;do perl reverse_complementary.pl -fa $i -out ${i%.*}-reverse.fasta; done

#然后重排（按照参考基因组）
for i in `ls *fasta`
do
python3 /work/users/wangyw/biosoft/release_MitoZ_v2.3/useful_scripts/Mitogenome_reorder.py \
-f $i \
-r /work/users/wangyw/02_S.micro/input/mt/S.micro.fa
done

----------------------------------------------------------------------
#还有一些需要换用小的kmer来组装
	#没有成环：
	JYP-03 --kmer51
	SCP-06 --kmer51
	SMHP-02 --kmer51
	ZZP-02.fasta --kmer51
	#注释有问题：
	FYP-03：notfound-trnT --kmer51
	JYP-15：Split/duplicated features: cob --kmer31

----------------------------------------------------------------------
#解压
for i in `ls`;do mkdir ../anno/${i%.*};done
for i in `ls`;do unzip $i -d ../anno/${i%.*};done
#把默认的文件名改为样本名
for i in `ls`;do cd $i; cd `ls`;cp result.fas ../../../anno-fas/$i.fas;cd ../..;done

## 为序列名填加上样本名
for i in `ls`;do cat $i|sed "s/>/>${i%.*};/g" > ../anno-fas2/$i ;done

----------------------------------------------------------------------

