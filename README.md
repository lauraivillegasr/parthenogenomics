### **ESTIMATION OF MUTATION RATES**



Use sanger encoding for all files (changing those with solexa): 

	cat rawreads.solfastq | seqret -filter -auto -sformat fastq-solexa -osformat fastq-sanger -out MALrawreads.fastq


Pre-processing 
1. Trimming

```fastp -i MALrawreadsforward.fq.gz  -I MALrawreadsreverse.fq.gz -o MALreadsforward.fa.gz -O MALreadsreverse.fa.gz -h reportfromreads```

2. Mapping

```bwamem2 mem -M -t 30 -R "@RG\tID:ASEX\tSM:PS1159_c12_PS_8\tPL:ILLUMINA\tPU:1" referencegenome.fa MALreadsreverse.fa.gz > MALmappedreadstoreference_bwamem.sam```

3. Convert sam to bam

```ls -1 | sed 's/_bwamem.sam//g' > list-XX```
```while read f; do samtools view -b $f"_bwamem.sam" > $f".bam" ;done < list-XX```

4. Sort files 

```samtools sort -o MALmappedreadstoreference.sort.bam MALmappedreadstoreference.bam```

4.1. Check coverage from sorted files

```ls *.sort.bam | parallel -j 5 'samtools depth {} > {}.sort.bam.depth'```

5. Remove duplicates using picard

```java -jar picard.jar MarkDuplicatesWithMateCigar I=MALmappedreadstoreference.sort.bam O=MALmappedreadstoreference.rmd.bam M=mappedreadstoreference.bam.metrics VALIDATION_STRINGENCY=SILENT MINIMUM_DISTANCE=300 REMOVE_DUPLICATES=true```

6. Remove low quality reads 

```ls *.sort.rmd.bam | parallel 'samtools view -q 30 -b {} > {.}.q30'```

Samtools flagstat can be used to check quality of mapping 

7. To check that the header is correctly stablished

```samtools view -H MALmappedreadstoreference.rmd.bam | grep '^@RG'```
```for f in *q30.bam ; do samtools view -H $f | grep '^@RG'; done```

8. Merging files for accuMUlate 

```samtools merge -r MALmergedmappedreads.sort.rmd.q30.bam MAL1mappedreadstoreference.sort.rmd.q30.bam MAL2mappedreadstoreference.sort.rmd.q30.bam MAL3mappedreadstoreference.sort.rmd.q30.bam MAL4mappedreadstoreference.sort.rmd.q30.bam MAL5mappedreadstoreference.sort.rmd.q30.bam MAL6mappedreadstoreference.sort.rmd.q30.bam```


9. Prepare data for accumulate, obtain ini file and obtain GC content using accumulate tools 


```samtools view -H MALmergedmappedreads.sort.rmd.q30.bam | python3 accuMulate-tools/extract_samples.py MAL_refpool - >> params.MAL.ini```
```python3 accuMulate-tools/GC_content.py referencegenome.fa >> params.MAL.ini```

	python3 dictionary_converter.py referencegneome.fa  > referencegenome.genomic.accu.dict


10. Generating windows using bedtools 

```mkdir -p /scratch/lvilleg1/MAL/Final_preprocessing/tmp```

	bedtools makewindows -g referencegenome.genomic.accu.dict -w 100000 | split -l 1 - tmp

11. Running accuMUlate to obtain candidate mutations 


```parallel -j 6 params.MAL.ini -b MALmergedmappedreads.sort.rmd.q30.bam -r referencegenome.fa -i {} ::: tmp/* > MAL_mutationcandidates``` 


12. Filtering according to different parameters to only keep mutations with high support

12.1. Define coverage range: 

With samtools depth we obtain a file with coverage at several positions, on R we can get the summary statistics for knowing the lower and upper range. We used 2 times the standard deviation of the ref pool for its upper limit: sd(file$V3)

The value for $11 changed according to the coverage range defined for the specific data set. 

12.2. Filter coverage range ($11), probability of having a mutation on a given site ($7, $8, $9), filter for unique mutations on a sample (not present in other lines $15), avoid calling a mutation given mismapped reads ($16 and $17), enough reads support the mutation ($18 and $19).

```cat MAL_mutationcandidates | awk '{if ($11 >=332 && $11 <=575 && $15 ==0 && $7 >=0.9 && $8 >=0.9 && $9 >=0.9 && $16 <=1.96 && $17 <=1.96 && $18 >=0.05 && $19 >=0.05) print $0}' > PS1159_mutationcandidates.filter-A.bed```

12.2.1. Step by step the decrease in number of putative mutations can be tracked. 

	cat PS1159_mutationcandidates | awk '{if ($11 >=332 && $11 <=575) print $0}' | wc -l 

	cat PS1159_mutationcandidates | awk '{if ($11 >=332 && $11 <=575 && $15 ==0) print $0}' | wc -l

	cat PS1159_mutationcandidates | awk '{if ($11 >=332 && $11 <=575 && $15 ==0 && $7 >=0.9 && $8 >=0.9 && $9 >=0.9) print $0}' | wc -l

	cat PS1159_mutationcandidates | awk '{if ($11 >=332 && $11 <=575 && $15 ==0 && $7 >=0.9 && $8 >=0.9 && $9 >=0.9 && $16 <=1.96 && $17 <=1.96) print $0}' | wc -l


13. Obtain the number of callable sites within the defined coverage region. 

We used the already obtained depth files from ```samtools depth filename.bam > filename.depth```

```cat MALmappedreadstoreference.sort.bam.depth | awk '{if ($3 >=10 && $3 <=50) print $0}' | wc -l ```


14. Obtaining mutation rates and confidence intervals: 

    15.1. For mutation rates the following equation was used for each of the mutation lines: 

	μ=(called mutations)/(generations∗callable sites)

    14.2. Average of μ for each of the strains was calculated. (Can’t insert equation, basically all μ divided the total number of μ for the strain). 

    14.3. Estimation of confidence intervals: 

Downloading Bayesian first aid on R

	install.packages("devtools")
	devtools::install_github("rasmusab/bayesian_first_aid")

	MAL_sex = c(1.62637E-09, 6.4472E-10, 1.60744E-09, 5.93941E-10, 1.17891E-09)
	MAL_asex =c(9.64506E-10, 6.94089E-10, 4.08937E-10, 7.2004E-10, 7.9846E-10)
	comparing= c(9.41896E-10, 5.67965E-10)
	sites = c(278665301,374228169)
	bayes.poisson.test(comparing, sites)
	plot(bayes.poisson.test(comparing, sites))

Sites refers to the number of callable sites for each reproduction mode. The result of this analysis tells us how different our groups are and provides confidence intervals of the values: a lower limit and an upper limit.



### **POPULATION ANALYSIS**

**Part 1.** **PRE-PROCESSING of data**
commands implemented on re-sequencing data from Illumina sequencing. 

1. Trimming using fastp: 

```fastp -i forward.fq.gz -I reverse.fq.gz -o out.forward.fastq -O out.reverse.fastq```


2. Mapping to reference genome via bwa mem

2.1 Create index for mapping: 
```bwamem2 index ref.fa```

2.2 Mapping: 
```bwamem2 mem -M -t 30 -R “@RG\tID:sample-id\tSM:sample\tPL:ILLUMINA\tPU:1” reference-genome.fasta out.forward.fastq out.reverse.fastq > popX_bwamem.sam``` 


3. Creating list for looping further analysis
  
```ls -1 | sed 's/_bwamem.sam//g' > list-XX```

4. Convert sam to bam and sort files


```ls -1 | sed 's/_bwamem.sam//g' > list-XX```
```while read f; do samtools view -b $f"_bwamem.sam" > $f".bam" ;done < list-XX```
```cat list-XX | parallel -j 12 'samtools sort -@ 4 -o {}bwamem.sort.bam {}.sam'```


5. Indexing sorted bam files 

 ```ls *.sort.bam | parallel samtools index '{}'```
	
6. Quality statistics of your mappings via qualimap. 
	#create data-description-file for multi-bamqc
```ls -d -1 $PWD/** | grep sort.bam$ > qualimap.list  ```

	#example:

|NAME        | PATH	|GROUP |  
-------------|----------|-------|
|pop1|	pop1_bwamem.sort.bam |      pop1|
|po2|	pop2_bwamem.sort.bam   |     pop2|
|popX|	popX_bwamem.sort.bam    |    popX|



	qualimap multi-bamqc -r -d qualimap.list --java-mem-size=2400M 

7. Removing duplicates with picard tools. 

		while read f; do java -jar picard.jar MarkDuplicatesWithMateCigar I=$f".sort.bam" O=$f".sort.rmd.bam" M=$f"sort.bam.metrics" VALIDATION_STRINGENCY=SILENT MINIMUM_DISTANCE=300 REMOVE_DUPLICATES=true & done < list


8. Filtering the data to keep just quality higher than 30:

```samtools view -q 30 -b popX_bwamem.sort.bam > popX_bwamem.sort.q30.bam```

9. Gatk haplotypecaller. 


    9.1. Indexing files after filtering
```ls *.bam | parallel samtools index '{}'```

    9.2. Haplotype calling 
```gatk HaplotypeCaller -R reference_genomes/panagrolaimus_es5.PRJEB32708.WBPS15.genomic.fa -I filesq30/P_bromber.sort.rmd.q30.bam -O P-bromber.vcf```


**ADDITIONAL INFORMATION ON PRE-PROCESSING:**

1. For p_superbus I used NGM - reads shorter than 70bp. 

	1.1 Create index for mapping:
```ngm -r reference-sequence.fa```

 	1.2 Mapping: 
```ngm -r reference-sequence.fa -1 forward_reads.fq -2 reverse_reads.fq -o results.sam ```


2. For files with insert size smaller than double read lenght we used pear. 

```pear -f forwardread.sanfastq.gz -r reverseread.sanfastq.gz -o pearoutputname```
				
pear output gives: file.unassembled.forward.fastq, file.unassembled.reverse.fastq, file.assembled.forward.fastq file.assembled.reverse.fastq and discarded reads

Mapping for this files was done twice: once for assembly and once for unassembled reads and then they were merged. 

```bwamem2 mem -M -t 30 -R "@RG\tID:ASEX\tSM:PS1806\tPL:ILLUMINA\tPU:1" reference-genome.fasta file.unassembled.forward.fastq file.unassembled.reverse.fastq > mapped_unassembledpear.sam```  #no extra indexing required, as the same index from previous mapping could be used

After re-mapping in either case, all steps from pre-processing were followed again. 

**Part 2.** **POPULATION ANALYSIS USING POPOOLATION:**


The sorted final bam files are used as initial input for the tool. 

1. Creating pileup file for (A) individual files (further estimation of pi, theta) and (B) merged files (for Fst estimation)

```samtools mpileup popX_bwamem.sort.bam  > popX.pileup```
```samtools mpileup -B -b list-samtoolpileup_Xsex > Xsexpop.mpileup```

2. Creating syncronized files for further estimations 

```java -ea -Xmx7g -jar popoolation2_1201/mpileup2sync.jar --input Xsexpop.mpileup --output Xsexpop_java.sync --fastq-type sanger —min-qual 30 --threads 4```

2.1. Fst estimated on sliding window non overlapping 

```perl popoolation2_1201/fst-sliding.pl --input Xsexpop_java.sync  --output Xsexpop_w1kb.fst  --suppress-noninformative --min-count 4 --min-coverage 10 --max-coverage 80 --min-covered-fraction 0,5 --window-size 1000 --step-size 1000 --pool-size 3000```
	
About coverage ranges: for asex populations 10 - 80 and for sex pops 10-56 (meaning 23 is the covergae per gene copy, ase population are triploid whereas sexual ones are diploid)
for MALpops a coverage range of 10 - 55 for asexual populations and 10-40 for the sexual populations were selected (15X per genome copy)

2.1.2 Extracting values of Fst from all columns - needs to be done for all columns with pairwise comparisson

	awk '{ gsub(/1:2=/,"", $6); print } ' Xsexpop_w1kb.fst > extractedpositioncolumnX

1:2 means fst of pop1 vs pop2. The order of populations is the same as provided on 1(B)

2.2. Fst estimation on single positions to obtain common positions between all populations (asex as a group and sex as another)

	perl popoolation2_1201/fst-sliding.pl --input Xsexpop_java.sync  --output Xsexpop_w1kperpos.fst  --suppress-noninformative --min-count 4 --min-coverage 10 --max-coverage 80 --min-covered-fraction 0,5 --window-size 1 --step-size 1 --pool-size 3000

2.2.1 From fst per portion (2.2), grab only the first two columns that have the information from positions present in all pops

	awk '{print $1, $2 }' Xsexpop_w1kperpos.fst > commonpositions_Xsex
			
3. With he positions file obtained in 2.2.1, common positions between all populations were extracted on the individual pileup files using: 

```awk 'NR==FNR{a[$1,$2]; next} ($1,$2) in a'  commonpositions_Xsex Xsexpop.mpileup > Xsexpop.corrected.mpileup```

3.1. Watterson theta and pi estimation using corrected pileup files

	perl popoolation_1.2.2/Variance-sliding.pl --measure theta --input Xsexpop.corrected.mpileup --output Xsexpop.mpileup_WT --pool-size 3000 --min-count 2 --min-coverage 10 --max-coverage 80 --window-size 1 --step-size 1 --fastq-type sanger
 
	perl popoolation_1.2.2/Variance-sliding.pl --measure pi --input Sex_network/p_davidi.corrected.pileup --output Sex_network/P_davidi_pi.file --pool-size 3000 --min-count 2 --min-coverage 10 --max-coverage 80 --window-size 1 --step-size 1 --fastq-type sanger

3.1.3 Remove all “fake/empty” positions, created by population, to only have meaningful results

	awk 'NR==FNR{a[$1,$2]; next} ($1,$2) in a' ommonpositions_Xsex Xsexpop.mpileup_WT > Xsexpop.mpileup_WT.corrected

3.1.4 Obtaining average values for watterson theta and pi 

	awk '{ total += $5; count++ } END { print total/count }' Xsexpop.mpileup_WT.corrected


Plots for visualizing results ere generated using ggplot2 on R.

4. Alternative theta estimation using Tetmer - This tool estimates theta on only homologous copies in (polyploid) our case triploid genomes, which allows us to see how the third copy porivdes more genetic diversity.

K-mer spectra of reads for each population was obtained as: 

	kat hist -o kmerspectra_fortetmer rawreadsforward.fq.gz rawreadsreverse.fq.gz
	
Provide histogram to tetmer and obtain per-k-mer theta, the autofit optio was selected for triploid allopoliploids. Since Kmer size k=27, the per Kmer theta result with the tool was divided by 27 to obtain the per nulceotide estimation. 

**Part 3.** **GENE NETWORK USING BUSCO GENES**

1. Run busco on the reference genomes to obtain complete genes that - done in Gvolante

| Characteristic| Description |
| --- | --- |
|Fasta file|preference_genome.fasta|
|Cut-off length for sequence statistics and composition|1|
|Sequence type|genome|
|Selected program |BUSCO_v4|
|Selected ortholog set |Nematoda|


2. With the list obtained in 1 - blast output coordinates.tsv as a list of genes, the regions of interest were extracted from each bam file using (

Editing coordinate file to concatenated name  (e.g. PSS159_contig1675 8 1943 —> PS1159_contig1675:8-1675) on Excel using command: ```(=B1&":"&C1&"-"&D1)```.
```while read f; do samtools view -b Xpop.sort.rmd.q30.bam $f >> $f".Xpop.bam" ; done <list_genes ```


3. Create pileup files with bcftools (summarising the base calls of aligned reads to a reference sequence). The -A refers to paired-end data, needs to be specified as it isn’t default. 

```while read f, do bcftools mpileup -A -f referencegenome.fa.gz $f > $f"_mpileup"; done <list_bamsort```

4. Variant calling using pileup files. -Ov option will give an uncompressed vcf file as result, other options are available such as BCF of vcf.gz, but vcf is required for following steps. 


```for n in *_mpileup ; do bcftools call -mv -Ov $n > ${n%*_mpileup}.vcf; done```

Note: usually steps 3 and 4 are performed as one command ```bcftools mpileup -Ou -f reference.fa alignments.bam | bcftools call -mv -Ob -o calls.bcf```

5. Create and index for each vcf file

```for file in *.vcf; do gatk-4.2.3.0/gatk IndexFeatureFile -I ${file} -O ${file}.idx ; done```

6. “Manual” editing. GATK wouldn’t take regions starting at position 0, but needed to be provided as position 1. All files were renamed when their position started with :0-.e.g. PS1159_contig1345:0-1928.ZZZ.ZZZ. 

```for f in *:0-*; do mv -i -- "$f" "${f//:0-/:1-}"; done ```

7. Create a list with all positions now with all positions starting with 1. Or change same “word” in list used in 2. 

8. Create index of fasta file 

```samtools faidx referencegenoem.genomic.fa```

9. Create a dictionary for the reference genome

```gatk-4.2.3.0/gatk CreateSequenceDictionary -R referencegenome.genomic.fa ```

10. Create a consensus sequence, it takes the reference genome we mapped against -R, provide the output name we want for our file -O, the region from the reference genome that we want to have a consensus sequence from while using the reads present in this region -L and the VCF file  

(the list only contains the contain position as it is common between all strains, individual scripts were submitted for each strain to allow manual specification of files to be used as input e.g. $j.bcf.vcf where region is $j

	while read f ;do gatk-4.2.3.0/gatk FastaAlternateReferenceMaker -R referencegenome.genomic.fa -O $f.Xpop.consensus.fa -L $f -V $f.Xpop.vcf ; 
	

11. Prepare files for the alignment, the header of each consensus sequence needs to be renamed so we can track from which strain it is coming from once we do the alignment. This command runs inside the folder of each strain 

```for FILE in *.fa;do awk '/^>/ {gsub(/.fa?$/,"",FILENAME);printf(">%s\n",FILENAME);next;} {print}' $FILE >>changed_${FILE}; done```

12. Get list of all consensus sequences present in all strains

    12.1 Obtain list of fasta file consensus from each folder/strain

```ls Xpop/changed_* | sed 's/Xpop\///g'> list_Xpop```

Xpop can be change in each case for the name of folder/strain

12.2 All lists where transferred to excel as one column, manual editing to just have the information on gene location was done using find&replace, then the duplicate function was used and only those positions present in ALL columns were kept and a new file was created: 

(e.g. 375_genes)

12.3 With a loop using the file created in 12.2, we copied all gene files from each strain to a new folder 

```while read f; do cp try/$f* try_2/$f.Xpop; done < 375_genes```

12.4 Put together all the files from one gene into one file for the mafft alignment

 ```while read f; do cat $f* >> $f"_for_mafft" ; done < 375_genes```

12.5 Check some files to make sure that you have all 7 headers from the 7 strains 

```ls *mafft | head -n 10```

```grep -o '>' Pops_contig10078:2915-14059_for_mafft | wc -l ```

You can also check the headers themselves 

```grep '>' changed_PS1159_contig10078:2915-14059_for_mafft```

The output looks like: 7 lines each with the name of the position the files refers to swell as the strain to which a sequence corresponds

	>PS1159_contig10078:2915-14059.DL137G2.consensus
	>PS1159_contig10078:2915-14059.JB051.consensus
	>PS1159_contig10078:2915-14059.p_davidi.consensus
	>PS1159_contig10078:2915-14059.PS1159.consensus
	>PS1159_contig10078:2915-14059.PS1162.consensus
	>PS1159_contig10078:2915-14059.PS1579.consensus
	>PS1159_contig10078:2915-14059.PS1806.consensus

All files where moved to a folder named “mafft_files”. 
 
13. Generate allignment

```parallel -j 8 'mafft {} > {.}.aligned.fasta' ::: *_for_mafft```

14. Concatenate the allignements and save the file in nexus format to make a split network

Using genious, the allignements per gene where uploaded and concatenated using

	Tools → Concatenate Sequences or Alignments
	
and exporting the resulting file in .nex format. 

15. Upload the nexus format file into SplitsTree and create a network using: ```Networks -> MedianNetwrok```
