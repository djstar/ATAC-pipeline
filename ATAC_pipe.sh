#!/bin/bash
#set -e
## ATAC pipeline 22 Aug 2017. Works on hg19.
# Make sure samtools is latest
source ~/.profile
#gunzip *.gz
while read p; do echo $p trim; trim_galore --paired $p\.R1.fastq.gz $p\.R2.fastq.gz;done < names.txt
# Map to hg19
while read p; do echo $p map hg19; bowtie2 -t -p 24 -N 1 --mp 2,1  --no-discordant --no-mixed --very-sensitive-local -I 100 -X 3000 -x ~/TopHatGenomes/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome -1 $p\.R1\_val\_1.fq.gz -2 $p\.R2\_val\_2.fq.gz 2>>mapping.stats| samtools view -@ 6 -b -S -o $p.paired.bt2.hg19.bam -;done < names.txt

# Remove blacklist reads
for name in *.bam; do echo $name Remove Blacklist; samtools sort -@ 8 -n $name | pairToBed -type neither -abam stdin -b ~/blacklist/wgEncodeDacMapabilityConsensusExcludable.bed > $(echo $name | sed 's/.bam//g').noBL.bam; rm $name; done;
# Sort reads
for name in *.bam; do echo $name; samtools sort -@ 16 $name -o $(echo $name | sed 's/.bam//g').sorted.bam; samtools index $(echo $name | sed 's/.bam//g').sorted.bam ; rm $name; done;

## Get ATAQC  done
mkdir ATACQC
cd ./ATACQC
ln -s ../*.bam* .
for name in *.bam; do echo $name Making ATAC pipeline; make_ataqc_pipeline $name; done;
# Run ATAQC pipeline
for name in qc*.bam.sh; do echo name; source ~/.profile; bash ./$name; done;
# Need this to make the html work
mkarv allATACQCout.web *.json --force

## Sort, index and make TDF
for name in *hqaa*.bam; do echo $name >> total_reads ; samtools view -F 8 $name | wc -l >> total_reads; done;
for name in *hqaa.bam; do echo $name INDEX; samtools index $name; done;
for name in *hqaa.bam; do echo $name Making TDF; ~/Programs/IGVTools/igvtools count -z 7 -w 25 $name $name.count.tdf hg19; done;
mkdir TDF
mv *.tdf ./TDF

# Shift peaks, convert from bam to bed
for name in *hqaa.bam; do echo $name shift cleavage sites; perl /data/home/sjiang/Programs/PIPES/bam2bed_shift.pl $name; done;

# Make bedgraph file and normalize
for name in *hqaa.bed; do echo $name Making Bedgraph; genomeCoverageBed -bg -split -i $name -g /data/home/sjiang/Programs/bedtools2/genomes/human.hg19.genome > $(echo $name | sed 's/.bed//g').bedGraph; rm $name; done;
for name in *hqaa.bedGraph; do echo $name Normalizing Bedgraph; python /data/home/sjiang/Programs/PIPES/normalize_bedGraph.py --bg $name --bam $(echo $name | sed 's/.bedGraph//g').bam > $(echo $name | sed 's/.bedGraph//g').norm.bedGraph; done;
# Make union pooled bed files (union regions from all ATAC peaks, All.ATAC.pooled.broadPeak)
# ATACQC are called *.hqaa.bam.macs2_peaks.broadPeak
#	1.	Make a combined broadpeak and take only the first 3 columns
cat *.hqaa.bam.macs2_peaks.broadPeak > ATAC.all.pooled.tmp.broadPeak
cut -f 1-3 ATAC.all.pooled.tmp.broadPeak > ATAC.all.pooled.broadPeak;

#	2.	Sort the bed files before merge
sortBed -i ATAC.all.pooled.broadPeak > ATAC.all.pooled.sorted.broadPeak
#	3.	Merge files, all peaks within 50bp of each other will be merged together
mergeBed -d 50 -i ATAC.all.pooled.sorted.broadPeak > All.ATAC.pooled.broadPeak
rm ATAC.all.pooled.*broadPeak
# Make intersecting bed files (intersecting regions common amounts ATAC peaks, or All.ATAC.common.broadPeak)
intersectBed -a All.ATAC.pooled.broadPeak -b *.hqaa.bam.macs2_peaks.broadPeak > All.ATAC.common.dup.broadPeak
sortBed -i All.ATAC.common.dup.broadPeak > All.ATAC.common.dup.sorted.broadPeak
mergeBed -i All.ATAC.common.dup.sorted.broadPeak > All.ATAC.common.broadPeak
rm *dup*.broadPeak
# Plot enrichment around bed files, and TSS, and pairwise correlations

# First, sort bedGraphs (only normalized ones!)
for name in *norm.bedGraph; do echo $name sorting BedGraph; sort -k1,1 -k2,2n $name > $(echo $name | sed 's/.bedGraph//g').sorted.bedGraph; rm $name; done;
# Convert BedGraphs to bigwigs
for name in *norm*bedGraph; do echo $name convert BedGraph to BigWig; bedGraphToBigWig $name ~/bowtieGenomes/genome_table.human.hg19.txt $(echo $name | sed 's/.bedGraph//g').bw; done

# Build matrix for union peaks
computeMatrix scale-regions -S *.bw -R All.ATAC.pooled.broadPeak --sortRegions descend --sortUsing mean --beforeRegionStartLength 5000 --regionBodyLength 1000 --afterRegionStartLength 5000 --blackListFileName ~/blacklist/wgEncodeDacMapabilityConsensusExcludable.bed -p 24 -o All.ATAC.pooled.regions_5-1-5kb.mat.gz

# Build matrix for intersect peaks
computeMatrix scale-regions -S *.bw -R All.ATAC.common.broadPeak --sortRegions descend --sortUsing mean --beforeRegionStartLength 5000 --regionBodyLength 1000 --afterRegionStartLength 5000 --blackListFileName ~/blacklist/wgEncodeDacMapabilityConsensusExcludable.bed -p 24 -o All.ATAC.common.regions_5-1-5kb.mat.gz

# Heatmap to visualize
plotHeatmap -m All.ATAC.pooled.regions_5-1-5kb.mat.gz --sortRegions descend --sortUsing mean --colorList "blue,white,red" --startLabel Start --endLabel End --plotTitle "Mean Descend Pooled peaks" -out All.ATAC.pooled.regions_5-1-5kb.Try1.png
plotHeatmap -m All.ATAC.common.regions_5-1-5kb.mat.gz --sortRegions descend --sortUsing mean --colorList "blue,white,red" --startLabel Start --endLabel End --plotTitle "Mean Descend Common peaks" -out All.ATAC.common.regions_5-1-5kb.Try1.png

# Do multibw summary
multiBigwigSummary BED-file -b *.bw -bl ~/blacklist/wgEncodeDacMapabilityConsensusExcludable.bed -p 24 -out multiBW_All.ATAC.pooled.regions.ngz --BED All.ATAC.pooled.broadPeak
multiBigwigSummary BED-file -b *.bw -bl ~/blacklist/wgEncodeDacMapabilityConsensusExcludable.bed -p 24 -out multiBW_All.ATAC.common.regions.ngz --BED All.ATAC.common.broadPeak

# Plot correlations between samples for 1. Heatmap 2. Scatterplot
plotCorrelation -in multiBW_All.ATAC.pooled.regions.ngz --corMethod pearson --skipZeros --plotTitle "Pearson Correlation between ATAC-MIBI pooled samples" --whatToPlot scatterplot -o scatterplot_pearson.All.ATAC.pooled.png
plotCorrelation -in multiBW_All.ATAC.pooled.regions.ngz --corMethod pearson --skipZeros --plotTitle "Pearson Correlation between ATAC-MIBI pooled samples" --whatToPlot heatmap -o heatmap_pearson.All.ATAC.pooled.png

plotCorrelation -in multiBW_All.ATAC.common.regions.ngz --corMethod pearson --skipZeros --plotTitle "Pearson Correlation between ATAC-MIBI common samples" --whatToPlot scatterplot -o scatterplot_pearson.All.ATAC.common.png
plotCorrelation -in multiBW_All.ATAC.common.regions.ngz --corMethod pearson --skipZeros --plotTitle "Pearson Correlation between ATAC-MIBI common samples" --whatToPlot heatmap -o heatmap_pearson.All.ATAC.common.png

# Plot PCA
plotPCA -in multiBW_All.ATAC.pooled.regions.ngz -o PCA-All.ATAC.pooled.png -T "PCA of read counts at all ATAC-MIBI pooled samples"
plotPCA -in multiBW_All.ATAC.common.regions.ngz -o PCA-All.ATAC.common.png -T "PCA of read counts at all ATAC-MIBI common samples"

# Get fragment sizes
for name in *.hqaa.bam; do echo Plotting Fragment Lengths; Rscript ~/Programs/PIPES/plot.bam.flen.R $name $(echo $name | sed 's/.bam//g').fraglen; done;

# NGS PLOT
mkdir ngsplot
cd ngsplot
for name in ../*hqaa.bam; do echo $name  -1 $(echo $name | sed 's/.noBL.sorted.hqaa.bam//g')>> config.tss; done;
for name in config*; do echo $name; ngs.plot.r -G hg19 -R tss -C $name -O $name.plot -P 16 -T ATAC-Signal -L 2000; done
