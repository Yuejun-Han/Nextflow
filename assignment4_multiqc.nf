
ch_fastq = Channel.fromFilePairs(params.fastq_pattern)
ch_star = Channel.fromFilePairs(params.fastq_pattern)
ch_staridx = Channel.value(file(params.star_ref))
ch_ref = Channel.value(file(params.ref_fa))
ch_ref_flat = Channel.value(file(params.ref_flat))
ch_ribosomal_intervals = Channel.value(file(params.ribosomal_intervals))



process fastqc{
  tag "FASTQC ${sample}"
  cpus 1
  memory '1 GB'
  time 1.hour
  module 'java-jdk/1.10.0_1'
  module 'fastqc/0.11.7'

  input:
    tuple val(sample), file(fq_files) from ch_fastq

  output:
    file("${sample}.R*.{html,zip}") into ch_fastqc_results

  script:
    """
    fastqc $fq_files
    """
}

process STAR_alignment{
  tag "STAR ${sample}"
  cpus 12
  memory '48 GB'
  time 24.hour
  module  'gcc/6.2.0'
  module 'STAR/2.6.1d'

  input:
   tuple val(sample), file(fq_files) from ch_star
   file(staridx) from ch_staridx 

  output:
   path("*.bam") into ch_sort_index
   path("*Log.final.out") into ch_star_log

  script:
    """
    STAR \
    --runThreadN 12 \
    --genomeDir $staridx \
    --readFilesIn $fq_files \
    --outSAMtype BAM Unsorted \
    --outSAMunmapped Within \
    --outFileNamePrefix $sample \
    --outSAMattrRGline ID:${sample}\tPU:${sample}\tSM:${sample}\tPL:ILLUMINA\tLB:${sample} \
    --outSAMattributes NH HI AS nM NM \
    --quantMode GeneCounts
    """
}

process sort_index{
  tag "sort_index ${bam}"
  cpus 6
  memory '6 GB'
  time 24.hour
  module 'gcc/6.2.0'
  module 'samtools/1.10'
  publishDir "hw4_results/"

  input:
    file(bam) from ch_sort_index

  output:
    path("${bam}.sorted.bam") into ch_sorted_bam
    path("${bam}.sorted.bam.bai")

  script:
    """
    samtools sort \
      -O bam \
      -o ${bam}.sorted.bam \
      -@ 6 \
      $bam

    samtools index ${bam}.sorted.bam 
    """
}


process post_alignment{
  tag "post_alignment {sorted_bam}"
  cpus 2
  memory '4 GB'
  time 24.hour
  module 'gcc/6.2.0'
  module 'java-jdk/1.8.0_92'
  module 'picard/2.18.29'
  module 'R/3.6.1'

  input:
    file(sorted_bam) from ch_sorted_bam
    file(ref_fa) from ch_ref
    file(ref_flat) from ch_ref_flat
    file(intervals) from ch_ribosomal_intervals


  output:
    file("${sorted_bam}.validation_metrics") into ch_valid_metrics
    file("${sorted_bam}.gc_bias_metrics") into ch_gcbias_metrics
    file("${sorted_bam}.gc_bias_chart.pdf")
    file("${sorted_bam}.gc_bias_summary")
    file("${sorted_bam}.alignment_summary") into ch_alignment_metrics
    file("${sorted_bam}.rna_seq_metrics") into ch_rnaseq_metrics


  script:
    """
    java -Xmx4g -jar \${PICARD} \
      ValidateSamFile \
      INPUT=${sorted_bam} \
      OUTPUT=${sorted_bam}.validation_metrics \
      REFERENCE_SEQUENCE=${ref_fa} \
      MODE=SUMMARY \
      VALIDATE_INDEX=true \
      INDEX_VALIDATION_STRINGENCY=EXHAUSTIVE

    java -Xmx4g -jar \${PICARD} \
      CollectGcBiasMetrics \
      INPUT=${sorted_bam} \
      OUTPUT=${sorted_bam}.gc_bias_metrics \
      REFERENCE_SEQUENCE=${ref_fa} \
      CHART_OUTPUT=${sorted_bam}.gc_bias_chart.pdf \
      SUMMARY_OUTPUT=${sorted_bam}.gc_bias_summary \
      IS_BISULFITE_SEQUENCED=false \
      METRIC_ACCUMULATION_LEVEL=ALL_READS \
      METRIC_ACCUMULATION_LEVEL=READ_GROUP \
      ALSO_IGNORE_DUPLICATES=true

    java -Xmx4g -jar \${PICARD}\
      CollectAlignmentSummaryMetrics \
      INPUT=${sorted_bam} \
      OUTPUT=${sorted_bam}.alignment_summary \
      REFERENCE_SEQUENCE=${ref_fa} \
      MAX_INSERT_SIZE=100000 \
      EXPECTED_PAIR_ORIENTATIONS=FR \
      IS_BISULFITE_SEQUENCED=false \
      METRIC_ACCUMULATION_LEVEL=ALL_READS \
      METRIC_ACCUMULATION_LEVEL=READ_GROUP

    java -jar -Xmx4g \${PICARD} \
      CollectRnaSeqMetrics \
      INPUT=${sorted_bam} \
      OUTPUT=${sorted_bam}.rna_seq_metrics \
      REF_FLAT=${ref_flat} \
      STRAND_SPECIFICITY=NONE \
      RIBOSOMAL_INTERVALS=${intervals} \
      METRIC_ACCUMULATION_LEVEL=ALL_READS \
      METRIC_ACCUMULATION_LEVEL=READ_GROUP \
      ASSUME_SORTED=true

    """
}


process multiqc{
  cpus 1
  memory '1 GB'
  time 24.hour
  module 'gcc/6.2.0'
  module 'multiqc/1.8.0'
  publishDir "hw4_results/"

  input:
    path("*_fastqc.zip") from ch_fastqc_results.collect()
    path("*Log.final.out") from ch_star_log.collect()
    path("*") from ch_valid_metrics.collect()
    path("*") from ch_gcbias_metrics.collect()
    path("*") from ch_alignment_metrics.collect()
    path("*") from ch_rnaseq_metrics.collect()

  output:
    file "*_report.html"
 
  script:
    """    
    multiqc . -f
    """
}

