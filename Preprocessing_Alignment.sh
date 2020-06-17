#Differential Expression analysis pipeline


# Noyau
cp -n R695_RH_C000VK1_3_1_H3TGKBBXY.IND12.fastq.gz R695_RH_C000VK1_3_2_H3TGKBBXY.IND12.fastq.gz R695_RH_C000VK4_3_1_H3TGKBBXY.IND2.fastq.gz R695_RH_C000VK4_3_2_H3TGKBBXY.IND2.fastq.gz R695_RH_C000VK7_4_1_H3TGKBBXY.IND16.fastq.gz \
R695_RH_C000VK7_4_2_H3TGKBBXY.IND16.fastq.gz R695_RH_C000VKA_4_1_H3TGKBBXY.IND19.fastq.gz R695_RH_C000VKA_4_2_H3TGKBBXY.IND19.fastq.gz R695_RH_C000VKD_5_1_H3TGKBBXY.IND4.fastq.gz R695_RH_C000VKD_5_2_H3TGKBBXY.IND4.fastq.gz \
R695_RH_C000VKG_5_1_H3TGKBBXY.IND12.fastq.gz R695_RH_C000VKG_5_2_H3TGKBBXY.IND12.fastq.gz R695_RH_C000VKJ_6_1_H3TGKBBXY.IND2.fastq.gz R695_RH_C000VKJ_6_2_H3TGKBBXY.IND2.fastq.gz R695_RH_C000VKM_6_1_H3TGKBBXY.IND16.fastq.gz \
R695_RH_C000VKM_6_2_H3TGKBBXY.IND16.fastq.gz R695_RH_C000VKP_7_1_H3TGKBBXY.IND12.fastq.gz R695_RH_C000VKP_7_2_H3TGKBBXY.IND12.fastq.gz R695_RH_C000VKS_3_1_H3W2MBBXY.IND12.fastq.gz R695_RH_C000VKS_3_2_H3W2MBBXY.IND12.fastq.gz \
R695_RH_C000VKV_3_1_H3W2MBBXY.IND19.fastq.gz R695_RH_C000VKV_3_2_H3W2MBBXY.IND19.fastq.gz R695_RH_C000VKY_3_1_H3W2MBBXY.IND7.fastq.gz R695_RH_C000VKY_3_2_H3W2MBBXY.IND7.fastq.gz R695_RH_C000VL1_7_1_H3TGKBBXY.IND19.fastq.gz \
R695_RH_C000VL1_7_2_H3TGKBBXY.IND19.fastq.gz R695_RH_C000VL4_8_1_H3TGKBBXY.IND4.fastq.gz R695_RH_C000VL4_8_2_H3TGKBBXY.IND4.fastq.gz R695_RH_C000VL7_8_1_H3TGKBBXY.IND6.fastq.gz R695_RH_C000VL7_8_2_H3TGKBBXY.IND6.fastq.gz \
R695_RH_C000VLA_8_1_H3TGKBBXY.IND19.fastq.gz R695_RH_C000VLA_8_2_H3TGKBBXY.IND19.fastq.gz R695_RH_C000VLD_1_1_H33CHBBXY.IND7.fastq.gz R695_RH_C000VLD_1_2_H33CHBBXY.IND7.fastq.gz R695_RH_C000VLG_1_1_H33CHBBXY.IND12.fastq.gz \
R695_RH_C000VLG_1_2_H33CHBBXY.IND12.fastq.gz R695_RH_C000VLJ_2_1_H33CHBBXY.IND2.fastq.gz R695_RH_C000VLJ_2_2_H33CHBBXY.IND2.fastq.gz R695_RH_C000VLM_2_1_H33CHBBXY.IND16.fastq.gz R695_RH_C000VLM_2_2_H33CHBBXY.IND16.fastq.gz \
R695_RH_C000VLP_3_1_H33CHBBXY.IND5.fastq.gz R695_RH_C000VLP_3_2_H33CHBBXY.IND5.fastq.gz R695_RH_C000VLS_3_1_H33CHBBXY.IND7.fastq.gz R695_RH_C000VLS_3_2_H33CHBBXY.IND7.fastq.gz R695_RH_C000VLV_4_1_H33CHBBXY.IND6.fastq.gz \
R695_RH_C000VLV_4_2_H33CHBBXY.IND6.fastq.gz R695_RH_C000VLY_4_1_H33CHBBXY.IND19.fastq.gz R695_RH_C000VLY_4_2_H33CHBBXY.IND19.fastq.gz R695_RH_C000VM1_5_1_H33CHBBXY.IND4.fastq.gz R695_RH_C000VM1_5_2_H33CHBBXY.IND4.fastq.gz \
R695_RH_C000VM4_5_1_H33CHBBXY.IND6.fastq.gz R695_RH_C000VM4_5_2_H33CHBBXY.IND6.fastq.gz R695_RH_C000VM7_6_1_H33CHBBXY.IND5.fastq.gz R695_RH_C000VM7_6_2_H33CHBBXY.IND5.fastq.gz R695_RH_C000VMA_6_1_H33CHBBXY.IND2.fastq.gz \
R695_RH_C000VMA_6_2_H33CHBBXY.IND2.fastq.gz R695_RH_C000VMD_6_1_H33CHBBXY.IND16.fastq.gz R695_RH_C000VMD_6_2_H33CHBBXY.IND16.fastq.gz R695_RH_C000VMG_7_1_H33CHBBXY.IND19.fastq.gz R695_RH_C000VMG_7_2_H33CHBBXY.IND19.fastq.gz \
R695_RH_C000VMJ_8_1_H33CHBBXY.IND4.fastq.gz R695_RH_C000VMJ_8_2_H33CHBBXY.IND4.fastq.gz R695_RH_C000VMM_8_1_H33CHBBXY.IND12.fastq.gz R695_RH_C000VMM_8_2_H33CHBBXY.IND12.fastq.gz R695_RH_C000VMP_1_1_H3W2MBBXY.IND2.fastq.gz \
R695_RH_C000VMP_1_2_H3W2MBBXY.IND2.fastq.gz R695_RH_C000VMS_1_1_H3W2MBBXY.IND16.fastq.gz R695_RH_C000VMS_1_2_H3W2MBBXY.IND16.fastq.gz R695_RH_C000VMV_2_1_H3W2MBBXY.IND5.fastq.gz R695_RH_C000VMV_2_2_H3W2MBBXY.IND5.fastq.gz \
R695_RH_C000VMY_2_1_H3W2MBBXY.IND7.fastq.gz R695_RH_C000VMY_2_2_H3W2MBBXY.IND7.fastq.gz R695_RH_C000VNA_4_1_H3W2MBBXY.IND16.fastq.gz R695_RH_C000VNA_4_2_H3W2MBBXY.IND16.fastq.gz R695_RH_C000VND_4_1_H3W2MBBXY.IND5.fastq.gz \
R695_RH_C000VND_4_2_H3W2MBBXY.IND5.fastq.gz R695_RH_C000VNG_5_1_H3W2MBBXY.IND2.fastq.gz R695_RH_C000VNG_5_2_H3W2MBBXY.IND2.fastq.gz R695_RH_C000VNJ_5_1_H3W2MBBXY.IND16.fastq.gz R695_RH_C000VNJ_5_2_H3W2MBBXY.IND16.fastq.gz \
R695_RH_C000VNM_6_1_H3W2MBBXY.IND19.fastq.gz R695_RH_C000VNM_6_2_H3W2MBBXY.IND19.fastq.gz R695_RH_C000VNP_6_1_H3W2MBBXY.IND4.fastq.gz R695_RH_C000VNP_6_2_H3W2MBBXY.IND4.fastq.gz R695_RH_C000VNS_7_1_H3W2MBBXY.IND12.fastq.gz \
R695_RH_C000VNS_7_2_H3W2MBBXY.IND12.fastq.gz R695_RH_C000VNV_7_1_H3W2MBBXY.IND19.fastq.gz R695_RH_C000VNV_7_2_H3W2MBBXY.IND19.fastq.gz R695_RH_C000VNY_7_1_H3W2MBBXY.IND7.fastq.gz R695_RH_C000VNY_7_2_H3W2MBBXY.IND7.fastq.gz \
R695_RH_C000VO1_8_1_H3W2MBBXY.IND16.fastq.gz R695_RH_C000VO1_8_2_H3W2MBBXY.IND16.fastq.gz R695_RH_C000VO4_8_1_H3W2MBBXY.IND5.fastq.gz R695_RH_C000VO4_8_2_H3W2MBBXY.IND5.fastq.gz R695_RH_C000VO7_1_1_H55NYBBXY.IND7.fastq.gz \
R695_RH_C000VO7_1_2_H55NYBBXY.IND7.fastq.gz R695_RH_C000VOD_2_1_H55NYBBXY.IND5.fastq.gz R695_RH_C000VOD_2_2_H55NYBBXY.IND5.fastq.gz R695_RH_C000VOG_2_1_H55NYBBXY.IND7.fastq.gz R695_RH_C000VOG_2_2_H55NYBBXY.IND7.fastq.gz /media/hdd/RTEL1RNA/Data/Noyau

# Chromatine
cp R695_RH_C000VK2_3_1_H3TGKBBXY.IND5.fastq.gz R695_RH_C000VK2_3_2_H3TGKBBXY.IND5.fastq.gz R695_RH_C000VK5_3_1_H3TGKBBXY.IND7.fastq.gz R695_RH_C000VK5_3_2_H3TGKBBXY.IND7.fastq.gz R695_RH_C000VK8_4_1_H3TGKBBXY.IND6.fastq.gz R695_RH_C000VK8_4_2_H3TGKBBXY.IND6.fastq.gz R695_RH_C000VKB_5_1_H3TGKBBXY.IND2.fastq.gz R695_RH_C000VKB_5_2_H3TGKBBXY.IND2.fastq.gz R695_RH_C000VKE_5_1_H3TGKBBXY.IND16.fastq.gz R695_RH_C000VKE_5_2_H3TGKBBXY.IND16.fastq.gz R695_RH_C000VKH_6_1_H3TGKBBXY.IND5.fastq.gz R695_RH_C000VKH_6_2_H3TGKBBXY.IND5.fastq.gz R695_RH_C000VKK_6_1_H3TGKBBXY.IND7.fastq.gz R695_RH_C000VKK_6_2_H3TGKBBXY.IND7.fastq.gz R695_RH_C000VL2_7_1_H3TGKBBXY.IND2.fastq.gz R695_RH_C000VL2_7_2_H3TGKBBXY.IND2.fastq.gz R695_RH_C000VL8_8_1_H3TGKBBXY.IND12.fastq.gz R695_RH_C000VL8_8_2_H3TGKBBXY.IND12.fastq.gz R695_RH_C000VLE_1_1_H33CHBBXY.IND16.fastq.gz R695_RH_C000VLE_1_2_H33CHBBXY.IND16.fastq.gz R695_RH_C000VLH_2_1_H33CHBBXY.IND5.fastq.gz R695_RH_C000VLH_2_2_H33CHBBXY.IND5.fastq.gz R695_RH_C000VLK_2_1_H33CHBBXY.IND7.fastq.gz R695_RH_C000VLK_2_2_H33CHBBXY.IND7.fastq.gz R695_RH_C000VLN_3_1_H33CHBBXY.IND6.fastq.gz R695_RH_C000VLN_3_2_H33CHBBXY.IND6.fastq.gz R695_RH_C000VLQ_3_1_H33CHBBXY.IND19.fastq.gz R695_RH_C000VLQ_3_2_H33CHBBXY.IND19.fastq.gz R695_RH_C000VLT_4_1_H33CHBBXY.IND4.fastq.gz R695_RH_C000VLT_4_2_H33CHBBXY.IND4.fastq.gz R695_RH_C000VLW_4_1_H33CHBBXY.IND12.fastq.gz R695_RH_C000VLW_4_2_H33CHBBXY.IND12.fastq.gz R695_RH_C000VLZ_5_1_H33CHBBXY.IND2.fastq.gz R695_RH_C000VLZ_5_2_H33CHBBXY.IND2.fastq.gz R695_RH_C000VMB_6_1_H33CHBBXY.IND7.fastq.gz R695_RH_C000VMB_6_2_H33CHBBXY.IND7.fastq.gz R695_RH_C000VME_7_1_H33CHBBXY.IND6.fastq.gz R695_RH_C000VME_7_2_H33CHBBXY.IND6.fastq.gz R695_RH_C000VMH_7_1_H33CHBBXY.IND2.fastq.gz R695_RH_C000VMH_7_2_H33CHBBXY.IND2.fastq.gz R695_RH_C000VMK_8_1_H33CHBBXY.IND16.fastq.gz R695_RH_C000VMK_8_2_H33CHBBXY.IND16.fastq.gz R695_RH_C000VMN_8_1_H33CHBBXY.IND5.fastq.gz R695_RH_C000VMN_8_2_H33CHBBXY.IND5.fastq.gz R695_RH_C000VMQ_1_1_H3W2MBBXY.IND7.fastq.gz R695_RH_C000VMQ_1_2_H3W2MBBXY.IND7.fastq.gz R695_RH_C000VMT_1_1_H3W2MBBXY.IND6.fastq.gz R695_RH_C000VMT_1_2_H3W2MBBXY.IND6.fastq.gz R695_RH_C000VMW_2_1_H3W2MBBXY.IND19.fastq.gz R695_RH_C000VMW_2_2_H3W2MBBXY.IND19.fastq.gz R695_RH_C000VMZ_2_1_H3W2MBBXY.IND4.fastq.gz R695_RH_C000VMZ_2_2_H3W2MBBXY.IND4.fastq.gz R695_RH_C000VNB_4_1_H3W2MBBXY.IND6.fastq.gz R695_RH_C000VNB_4_2_H3W2MBBXY.IND6.fastq.gz R695_RH_C000VNH_5_1_H3W2MBBXY.IND7.fastq.gz R695_RH_C000VNH_5_2_H3W2MBBXY.IND7.fastq.gz R695_RH_C000VNK_5_1_H3W2MBBXY.IND6.fastq.gz R695_RH_C000VNK_5_2_H3W2MBBXY.IND6.fastq.gz R695_RH_C000VNN_6_1_H3W2MBBXY.IND2.fastq.gz R695_RH_C000VNN_6_2_H3W2MBBXY.IND2.fastq.gz R695_RH_C000VNQ_6_1_H3W2MBBXY.IND16.fastq.gz R695_RH_C000VNQ_6_2_H3W2MBBXY.IND16.fastq.gz R695_RH_C000VO2_8_1_H3W2MBBXY.IND6.fastq.gz R695_RH_C000VO2_8_2_H3W2MBBXY.IND6.fastq.gz R695_RH_C000VO5_8_1_H3W2MBBXY.IND19.fastq.gz R695_RH_C000VO5_8_2_H3W2MBBXY.IND19.fastq.gz R695_RH_C000VO8_1_1_H55NYBBXY.IND4.fastq.gz R695_RH_C000VO8_1_2_H55NYBBXY.IND4.fastq.gz R695_RH_C000VOB_1_1_H55NYBBXY.IND6.fastq.gz R695_RH_C000VOB_1_2_H55NYBBXY.IND6.fastq.gz R695_RH_C000VOE_2_1_H55NYBBXY.IND19.fastq.gz R695_RH_C000VOE_2_2_H55NYBBXY.IND19.fastq.gz R695_RH_C000VOH_2_1_H55NYBBXY.IND4.fastq.gz R695_RH_C000VOH_2_2_H55NYBBXY.IND4.fastq.gz /media/hdd/RTEL1RNA/Data/Chromatine

# Preprocessing

for f in *.fastq.gz
./Tools/FastQC/fastqc $f -o /home/ismail/RTEL1RNA/QC_report_fastq

# Rename:
for file in *; 
do
mv $file `echo $file | cut -c1-19`
done

for f in *; 
do
mv $f "$f.fastq.gz"
done

#Clipping 

java -jar Tools/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 70 -phred33 /media/hdd/RTEL1RNA/Data/Noyau/R695_RH_C000VK1_3_1.fastq.gz /media/hdd/RTEL1RNA/Data/Noyau/R695_RH_C000VK1_3_2.fastq.gz /media/hdd/RTEL1RNA/Data/Noyau/Paired/Preprocessed_data/R695_RH_C000VK1_3_1_paired.fastq.gz /media/hdd/RTEL1RNA/Data/Noyau/Unpaired/Preprocessed_data/R695_RH_C000VK1_3_1_unpaired.fastq /media/hdd/RTEL1RNA/Data/Noyau/Paired/Preprocessed_data/R695_RH_C000VK1_3_2_paired.fastq /media/hdd/RTEL1RNA/Data/Noyau/Preprocessed_data/Unpaired/R695_RH_C000VK1_3_2_unpaired.fastq.gz ILLUMINACLIP:/home/ismail/Tools/Trimmomatic-0.39/adapters/TruSeq2-PE.fa:2:30:10 MINLEN:101

#Loop
cd  /media/hdd/RTEL1RNA/Data/Chromatine/

for R1 in *_1.fastq.gz
do
	R2=${R1%%_1.fastq.gz}"_2.fastq.gz"
	R1Paired=${R1%%.fastq.gz}"_paired.fastq.gz"
	R1Unpaired=${R1%%.fastq.gz}"_unpaired.fastq.gz"
	R2Paired=${R2%%.fastq.gz}"_paired.fastq.gz"
	R2Unpaired=${R2%%.fastq.gz}"_unpaired.fastq.gz"

java -jar /home/ismail/Tools/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 70 -phred33 $R1 $R2 /media/hdd/RTEL1RNA/Data/Chromatine/Preprocessed_data/Paired/$R1Paired /media/hdd/RTEL1RNA/Data/Chromatine/Preprocessed_data/Unpaired/$R1Unpaired /media/hdd/RTEL1RNA/Data/Chromatine/Preprocessed_data/Paired/$R2Paired /media/hdd/RTEL1RNA/Data/Chromatine/Preprocessed_data/Unpaired/$R2Unpaired ILLUMINACLIP:/home/ismail/Tools/Trimmomatic-0.39/adapters/TruSeq2-PE.fa:2:30:10 MINLEN:101

done

#Alignment with hisat2

#Create the reference index

./Tools/hisat2-2.1.0/hisat2-build /home/ismail/RTEL1RNA/Data/Ref/hg38.fa RNA_RTEL1/Data/Ref/hg38

#loop

cd /media/hdd/RTEL1RNA/Data/Chromatine/Preprocessed_data/Paired/

for R1 in *_1_paired.fastq.gz
do
	R2=${R1%%_1_paired.fastq.gz}"_2_paired.fastq.gz"
	bam=${R1%%_1_paired.fastq.gz}".bam"
	
/home/ismail/Tools/hisat2-2.1.0/hisat2 -p 30 -x /home/ismail/RTEL1RNA/Data/Ref/hg38 -1 $R1 -2 $R2 | samtools view -Sb - > /media/hdd/RTEL1RNA/Data/Chromatine/Alignment/$bam
done


#Sort Bam
#loop

cd /media/hdd/RTEL1RNA/Data/Chromatine/Alignment/

for b in *.bam
do
	sortedbam=${b%%.bam}".sorted.bam"

java -Xmx8G -jar /home/ismail/Tools/picard.jar SortSam INPUT=$b OUTPUT=Sorted/$sortedbam SORT_ORDER=coordinate
done



#create the bam index

cd /media/hdd/RTEL1RNA/Data/Chromatine/Alignment/Sorted

for b in *.bam
do
	index=${b%%.bam}".bai"
samtools index $b $index
done

#Count matrix using FeatureCount
# If the input contains location-sorted paired-end reads, featureCounts will automatically re-order the reads to place next to each other the reads from the same pair before counting them.
# Sometimes name-sorted paired-end input reads are not compatible with featureCounts (due to for example reporting of multi-mapping results) and in this case featureCounts will also automatically re-order them

# Noyau

featureCounts -T 10 -p -a /home/ismail/RTEL1RNA/Data/Annotations/Homo_sapiens.GRCh38.99.gtf -o /home/ismail/RTEL1RNA/Data/Noyau/Count_matrices/New/count_matrix_exon_raw.txt /media/ismail/Data_analysis/Noyau/Alignment/*.bam 2> /home/ismail/RTEL1RNA/Data/Noyau/Count_matrices/New/Exon.txt

featureCounts -T 10 -p  -t gene  -a /home/ismail/RTEL1RNA/Data/Annotations/Homo_sapiens.GRCh38.99.gtf -o /home/ismail/RTEL1RNA/Data/Noyau/Count_matrices/New/count_matrix_genes_raw.txt /media/ismail/Data_analysis/Noyau/Alignment/*.bam 2> /home/ismail/RTEL1RNA/Data/Noyau/Count_matrices/New/Gene.txt

# Cytoplasme
featureCounts -T 10 -p -a /home/ismail/RTEL1RNA/Data/Annotations/Homo_sapiens.GRCh38.99.gtf -o /home/ismail/RTEL1RNA/Data/Cytoplasme/Count_matrices/New/count_matrix_exon_raw.txt /media/ismail/Data_analysis/Cytoplasme/Alignment/*.bam 2> /home/ismail/RTEL1RNA/Data/Cytoplasme/Count_matrices/New/Exon.txt

featureCounts -T 10 -p  -t gene  -a /home/ismail/RTEL1RNA/Data/Annotations/Homo_sapiens.GRCh38.99.gtf -o /home/ismail/RTEL1RNA/Data/Cytoplasme/Count_matrices/New/count_matrix_genes_raw.txt /media/ismail/Data_analysis/Cytoplasme/Alignment/*.bam 2> /home/ismail/RTEL1RNA/Data/Cytoplasme/Count_matrices/New/Gene.txt

#Chromatine
featureCounts -T 10 -p  -t gene  -a /home/ismail/RTEL1RNA/Data/Annotations/Homo_sapiens.GRCh38.99.gtf -o /home/ismail/RTEL1RNA/Data/Chromatine/Count_matrices/count_matrix_genes_raw.txt /media/hdd/RTEL1RNA/Data/Chromatine/Alignment/*.bam 2> /home/ismail/RTEL1RNA/Data/Chromatine/Count_matrices/summary.txt

