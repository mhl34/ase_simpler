#!/usr/bin/env bash
#SBATCH -e run_pipelineB136_T1_HS_error
#SBATCH -J run_pipelineB136_T1_HS.slurm
#SBATCH -o run_pipelineB136_T1_HS_output
#SBATCH --mail-type=END
#SBATCH --mail-user=mhl34@duke.edu
# --mem=120G
#SBATCH --mem=200G
#SBATCH --cpus-per-task=8

# modules
module load cutadapt/2.3-gcb01
module load picard-tools
module load java
module load samtools/1.3.1-gcb01
module load STAR
module load perl
module load Anaconda/4.3.0-fasrc01
module load python

# activate conda environment
source activate aseSimpler

# # directories
sampleName=B136_T1_HS
homeDir="/data/allenlab/sz_atac/data.transfer/to_Yuncheng"
directory="/data/allenlab/sz_atac/data.transfer/to_Yuncheng/${sampleName}"
supDirectory="/data/allenlab/sz_atac/data.transfer/to_Yuncheng/supFiles/"

# # echo "0) transfer file to homeDir"
# # DO NOT COMMENT OUT
# # scp -r mhl34@hardac-xfer.genome.duke.edu:/nfs/DHTS_crawfordlab/crawfordlab/rawData/human/psychENCODE_frontalcortex/ATACseq/${sampleName} /data/allenlab/sz_atac/data.transfer/to_Yuncheng

# echo "1) run_cutadapt.slurm"
# for FILE in ${directory}/*.fastq.gz; do
#     # cutadapt --quality-base 33 --minimum-length 17 -q 30 -a CTGTCTCTTATACACATCT -o ${directory}/*.trim.fastq ${directory}/*.fastq.gz
#     # input: ${directory}/${base}.fastq.gz
#     # output: ${directory}/${base}.trim.fastq
#     # function: cutadapt all of the fastq files
#     name=${FILE%.fastq.gz}
#     base=${name##*/}
#     echo "${base}.trim.fastq.gz"
#     cutadapt --quality-base 33 --minimum-length 17 -q 30 -a CTGTCTCTTATACACATCT -o "${directory}/${base}.trim.fastq.gz" "${directory}/${base}.fastq.gz"
# done

# echo "2) run_gzip.slurm"
# # for FILE in ${directory}/*.trim.fastq; do
# #     # input: ${directory}/${base}.trim.fastq
# #     # output: ${directory}/${base}.trim.fastq.gz
# #     # function: gzip all of the trimmed reads
# #    name=${FILE%.fastq.gz}
# #    base=${name##*/}
# #    echo "${base}.trim.fastq.gz"
# #    gzip "${directory}/${base}.trim.fastq"
# # done

# ulimit -n 100000

# echo "3) run_STAR.slurm"
# for FILE in ${directory}/*.trim.fastq.gz; do
#     # input: ${directory}/${base}.trim.fastq.gz
#     # output: ${directory}/${base}Aligned.sortedByCoord.out.bam
#     # function: STAR align all of the trimmed files onto hg19
#     name=${FILE%.trim.fastq.gz}
#     base=${name##*/}
#     echo "${directory}/${base}.trim.fastq.gz"
#     STAR --twopassMode Basic --runThreadN 24 --genomeDir /data/allenlab/sz_atac/data.transfer/to_Yuncheng/hg19_STAR_index         --readFilesIn "${directory}/${base}.trim.fastq.gz"         --alignEndsType EndToEnd         --waspOutputMode SAMtag         --varVCFfile "${supDirectory}/final288.all_chrm.vcf"         --outFilterMismatchNmax 10         --outSAMtype BAM SortedByCoordinate         --outReadsUnmapped Fastx         --outSAMattributes NH HI NM MD AS nM jM jI XS vA vG vW         --readFilesCommand "gunzip -c"         --outFileNamePrefix "${base}"
#     mv "${base}Aligned.sortedByCoord.out.bam" ${directory}
#     rm -rf ${base}*
# done
    
# echo "4) run_removeUnmapped.slurm"
# rm -rf "${directory}/bamList.txt"
# echo "#!/bin/bash" >> "${directory}/bamList.txt"
# printf '
# ' >> "${directory}/bamList.txt"
# echo -n "samtools merge ${directory}/${sampleName}merged.se.wq.bam " >> "${directory}/bamList.txt"

# for FILE in ${directory}/*Aligned.sortedByCoord.out.bam; do
#     # input: ${directory}/${base}Aligned.sortedByCoord.out.bam
#     # output: ${directory}/${base}.se.wq.bam
#     # function: remove the unmapped reads
#     name=${FILE%Aligned.sortedByCoord.out.bam}
#     base=${name##*/}
#     echo "${FILE}"
#     echo -n "${directory}/${base}.se.wq.bam " >> "${directory}/bamList.txt"
#     samtools view -b -F 1804 -q 30 "${FILE}" | samtools sort -@ 8 -O bam - > "${directory}/${base}.se.wq.bam"
# done

# echo "5) run_mergeBam.slurm"
# # input: ${directory}/${name}.se.wq.bam
# # output: ${directory}/${name}merged.se.wq.sorted.bam
# # function: merge the bam files
# name=${sampleName}
# samtools merge -b "${directory}/bamList.txt" -o "${directory}/${name}merged.se.wq.bam"
# chmod 701 "${directory}/bamList.txt"
# srun "${directory}/bamList.txt"
# samtools merge ${directory}/${name}merged.se.wq.bam ${directory}/*.se.wq.bam
# echo "${directory}/${name}merged.se.wq.bam"
# samtools sort -@ 8 -O bam "${directory}/${name}merged.se.wq.bam" > "${directory}/${name}merged.se.wq.sorted.bam"

# echo "6) run_dedup.slurm"
# # input: ${directory}/${name}merged.se.wq.sorted.bam
# # output: ${directory}/${name}merged.se.wq.sorted.dedup.bam
# # function: dedup the file
# java -jar $PICARD_TOOLS_HOME/picard.jar MarkDuplicates ASSUME_SORTED=true REMOVE_DUPLICATES=true I="${directory}/${name}merged.se.wq.sorted.bam"  O="${directory}/${name}merged.se.wq.sorted.dedup.bam" M=metrics.txt
# samtools index -b "${directory}/${name}merged.se.wq.sorted.dedup.bam"
# # mkdir "${homeDir}/alignment/${sampleName}"
# # cp "${directory}/${name}merged.se.wq.sorted.dedup.bam" "${homeDir}/alignment/${sampleName}"
# # cp "${directory}/${name}merged.se.wq.sorted.dedup.bam.bai" "${homeDir}/alignment/${sampleName}"

# echo "7) run_pileup.slurm"
# # input: ${directory}/${name}merged.se.wq.sorted.dedup.bam
# # output: ${directory}/${name}merged.se.final.sorted.pileup
# # function: pileup
# samtools faidx /data/allenlab/sz_atac/data.transfer/to_Yuncheng/hg19/hg19.fa
# het_sites_for_mpileup="${supDirectory}/final288.all_chrm.final.bed"
# ref="/data/allenlab/sz_atac/data.transfer/to_Yuncheng/hg19/hg19.fa"
# samtools mpileup -d 0 -B -s -f $ref --positions $het_sites_for_mpileup "${directory}/${name}merged.se.wq.sorted.dedup.bam" > "${directory}/${name}merged.se.final.sorted.pileup"

# # split the files into different chromosomes
# echo "8) run_pileup2base.slurm"
# # input: ${directory}/${name}.se.final.sorted.pileup
# # output: ${directory}/${name}output.txt
# # function: get by base reads
# perl pileup2base.pl "${directory}/${name}merged.se.final.sorted.pileup" 0 "${directory}/${name}output.txt"


echo "9) split files"
readCounts="${directory}/${sampleName}output.txt"
allChrmVCF="${supDirectory}/final288.all_chrm.vcf_abbrev.txt"
peaks="${supDirectory}/quant_peaks_1FDR_minOverlap2_300bpExt_150_TMM_no_rep_outliers.mid_peak.288.txt"
outDir="${directory}"
mkdir "${outDir}/output/"
mkdir "${outDir}/vcf/"
mkdir "${outDir}/peaks/"
python splitFileByChrom.py -rc "${readCounts}" -v "${allChrmVCF}" -p "${peaks}" -o "${outDir}"

echo "10) run_aseSimpler.slurm"
# input: readCounts, allChrmVCF, peaks
# output: p values
# function: get by base reads
for FILE in ${directory}/vcf/*;
do
    file="${FILE}"
    filename="${file##*/}"
    filename="${filename%.*}"
    readCounts="${directory}/output/${filename}.txt"
    allChrmVCF="${directory}/vcf/${filename}.txt"
    peaks="${directory}/peaks/${filename}.txt"
    outDir="${directory}"
    python ase_simplerFinal.py -rc "${readCounts}" -v "${allChrmVCF}" -p "${peaks}" -o "${outDir}" -ch "${filename}"
done

# cleanup
rm -rf ${directory}/*.fastq.gz
rm -rf ${directory}/*sortedByCoord.out.bam
rm -rf ${directory}/*se.wq.bam
rm -rf ${directory}/*merged.se.wq.sorted.bam
rm -rf ${directory}/output/
rm -rf ${directory}/peaks/
rm -rf ${directory}/vcf/
