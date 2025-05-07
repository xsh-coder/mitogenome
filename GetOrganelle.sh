#🌼🌼🌼extractfq 提取部分片段的数据(#提取的时候一般提取2个Gb即可，-size_required后面跟提取的数量)🌼🌼🌼
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
-----------------------------------------
###🌼🌼🌼利用提取的数据组装线粒体🌼🌼🌼
source activate geto
cd /data01/zhaoyh/FangYT/xush-mito
for i in *_1_extractfq.fq.gz; do
    filename=$(basename "$i")              # 取出文件名部分
    sample=${filename%_1_extractfq.fq.gz}  # 去掉后缀变成样本名
    get_organelle_from_reads.py \
        -1 ${sample}_1_extractfq.fq.gz \
        -2 ${sample}_2_extractfq.fq.gz \
        -s Sang-hap1.chr.fasta \
        -R 10 \
        -k 21,45,65,85,105 \
        -F animal_mt \
        -t 30 \
        -o ${sample}_out
done
#其中-1为5'端向3'端测得的序列；-2为3'端向5'端测得的序列，-o为输出的文件夹，-R为迭代循环次数，-k为打断成的kmer大小（为奇数，不能超过测序长度，例如pe150测序不能超过149），-f为组装的类型。
#另外，可以加入-s参数指定参考线粒体基因组（可以放入近缘种的做参考），-t参数来指定使用的cpu数量，更多参数可以去官网看一下。
#组装完成后结果序列文件中如果有circular即是成环。
#基本上这样就可以了，如果没成环可以调整一下-R、-k参数，或者用NOVOplasty以最长的contig为seed进行延伸（不建议用这个方法，因为可能有不定碱基或者gap的出现）。
#除非控制区特别长，或者样本有问题，一般二代就可以完成动物线粒体基因组组装。
-------------------------------------
###优先使用mitofish进行注释
#网页：https://mitofish.aori.u-tokyo.ac.jp/。选择中间的annotate(https://mitofish.aori.u-tokyo.ac.jp/annotation/input/),在下面选择需要导入的文件
#可以将后面的是否成环和可视化取消勾选，然后"Annotate"，下面有一个View Results，点击后有"Download MitoAnnotator Results"
#即可下载结果

