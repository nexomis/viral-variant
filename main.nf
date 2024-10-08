#!/usr/bin/env nextflow

nextflow.preview.output = true
include { validateParameters; paramsHelp; paramsSummaryLog; fromSamplesheet } from 'plugin/nf-validation'
include { showSchemaHelp; extractType } from './modules/config/schema_helper.nf'

log.info """
    |            #################################################
    |            #    _  _                             _         #
    |            #   | \\| |  ___  __ __  ___   _ __   (_)  __    #
    |            #   | .` | / -_) \\ \\ / / _ \\ | '  \\  | | (_-<   #
    |            #   |_|\\_| \\___| /_\\_\\ \\___/ |_|_|_| |_| /__/   #
    |            #                                               #
    |            #################################################
    |
    | viral-variant: Calling variants with peptide impact and comparison between samples (via transfer annotation): From raw reads, ref genome (with its own annotation or with other sequence and annotation to transfert annotation about PSA) and mapping of reads.
    |
    |""".stripMargin()

if (params.help) {
  log.info paramsHelp("nextflow run nexomis/viral-variant --input </path/to/samplesheet> [args]")
//  log.info showSchemaHelp("assets/input_schema.json")
//  log.info showSchemaHelp("assets/ref_genomes_schema.json")
  exit 0
}
validateParameters()
log.info paramsSummaryLog(workflow)

file(params.out_dir + "/nextflow").mkdirs()
// groovy fonction within nextflow script
def parse_sample_entry(it) {
  def type = "SR"
  def files = [file(it[1])]
  if (it[2] && !it[2].isEmpty() ) {
    files << file(it[2])
    type = "PE"
  } else {
    if (it[1].toString().toLowerCase().endsWith("spring")) {
      type = "spring"
    }
  }
  meta = [
    "id": it[0],
    "read_type": type,
    "do_primary": it[3],
    "batch_id": it[4],
    "rank_in_batch": it[5],
    "ref_id": it[6],
    "annot_id": (it[7] && !it[7].isEmpty() ) ? it[7] : null
  ]
  return [meta, files]
}
  
// include
include {PRIMARY_FROM_READS} from './modules/subworkflows/primary/from_reads/main.nf'
include {VIRAL_VARIANT} from './modules/subworkflows/viral_variant/main.nf'



workflow {
  inReads = Channel.of(
    [ [id: 'sars_n4_PeAn', batch_id: 'clinic', rank_in_batch: 5], [file('data/test/private_inputs/reads/sars_n4_PeAn_R1_subset.fastq.gz'), file('data/test/private_inputs/reads/sars_n4_PeAn_R2_subset.fastq.gz')] ],
    [ [id: 'sars_n6_XBB15', batch_id: 'batch1', rank_in_batch: 4], [file('data/test/private_inputs/reads/sars_n6_XBB15_R1_subset.fastq.gz'), file('data/test/private_inputs/reads/sars_n6_XBB15_R2_subset.fastq.gz')] ],
    [ [id: 'sars_n2_FrLa', batch_id: 'clinic', rank_in_batch: 0], [file('data/test/private_inputs/reads/sars_n2_FrLa_R1_subset.fastq.gz'), file('data/test/private_inputs/reads/sars_n2_FrLa_R2_subset.fastq.gz')] ]
  )
  //  [ [id: 'sars_n3_KeHe', batch_id: 'clinic', rank_in_batch: 2],['/mnt/d/viral_variant_test_dataset/inReads/sars_n3_KeHe_R1.fastq.gz', '/mnt/d/viral_variant_test_dataset/inReads/sars_n3_KeHe_R2.fastq.gz'] ]
  //)

  inRef = Channel.of(
    [ [id: 'batch1'], [file('data/test/private_inputs/ref/GCF_009858895.2_ASM985889v3_genomic.fna'), file('data/test/private_inputs/ref/GCF_009858895.2_ASM985889v3_genomic.gff')] ],
    [ [id: 'clinic'], [file('data/test/private_inputs/ref/sars_assembled.fasta')] ]
  )
  //  [ [id: 'batch1'], ['https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/858/895/GCF_009858895.2_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.fna.gz', 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/858/895/GCF_009858895.2_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.gff.gz'] ],

  inAnnot = Channel.of(
    [ [id: 'clinic'], [file('data/test/private_inputs/ref/GCF_009858895.2_ASM985889v3_genomic.fna'), file('data/test/private_inputs/ref/GCF_009858895.2_ASM985889v3_genomic.gff')] ]
  )
  //  [ [id: 'clinic'], ['https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/858/895/GCF_009858895.2_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.fna.gz', 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/858/895/GCF_009858895.2_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.gff.gz'] ]

  VIRAL_VARIANT(inReads, inRef, inAnnot)

  publish:
  VIRAL_VARIANT.out.var_by_batch >> 'var_by_batch'  
  VIRAL_VARIANT.out.var_by_smpl_raw >> 'var_by_smpl_raw'                        //
  VIRAL_VARIANT.out.subset_cov >> 'subset_cov'                                  //
  VIRAL_VARIANT.out.var_by_smpl_filter >> 'var_by_smpl_filter'                  //
  VIRAL_VARIANT.out.var_by_batch_ind_smpl_file >> 'var_by_batch_ind_smpl_file'
  VIRAL_VARIANT.out.transfered_gff >> 'transfered_gff'
  VIRAL_VARIANT.out.psa_algn >> 'psa_algn'                                      //
  VIRAL_VARIANT.out.psa_genomic_coords >> 'psa_genomic_coords'                  //
  VIRAL_VARIANT.out.flagstat >> 'flagstat'
  // bam in option: save_bam ?

}


output {
  directory 'out_dir'
  mode 'rellink'
}

