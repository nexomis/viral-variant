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
    [ [id: 'b1_P0', batch_id: 'batch1', rank_in_batch: 1], [file('data/test/input/reads/batch1/P0_R1.fq.gz'), file('data/test/input/reads/batch1/P0_R2.fq.gz')] ],
    [ [id: 'b1_P1', batch_id: 'batch1', rank_in_batch: 2], [file('data/test/input/reads/batch1/P1_R1.fq.gz'), file('data/test/input/reads/batch1/P1_R2.fq.gz')] ],
    [ [id: 'b2_P0', batch_id: 'batch2', rank_in_batch: 1], [file('data/test/input/reads/batch2/P0_R1.fq.gz'), file('data/test/input/reads/batch2/P0_R2.fq.gz')] ],
    [ [id: 'b2_P1', batch_id: 'batch2', rank_in_batch: 2], [file('data/test/input/reads/batch2/P1_R1.fq.gz'), file('data/test/input/reads/batch2/P1_R2.fq.gz')] ],
    [ [id: 'b2_P2', batch_id: 'batch2', rank_in_batch: 3], [file('data/test/input/reads/batch2/P1_R1.fq.gz'), file('data/test/input/reads/batch2/P1_R2.fq.gz')] ],
    [ [id: 'b2_P3', batch_id: 'batch2', rank_in_batch: 4], [file('data/test/input/reads/batch2/P1_R1.fq.gz'), file('data/test/input/reads/batch2/P1_R2.fq.gz')] ]
  )

  inRef = Channel.of(
    [ [id: 'batch1'], [file('data/test/input/ref/random.fa')] ],
    [ [id: 'batch2'], [file('data/test/input/ref/refseq.fa'), file('data/test/input/ref/refseq.gff')] ]
  )

  inAnnot = Channel.of(
    [ [id: 'batch1'], [file('data/test/input/ref/refseq.fa'), file('data/test/input/ref/refseq.gff')] ]
  )


  VIRAL_VARIANT(inReads, inRef, inAnnot)

  publish:
  VIRAL_VARIANT.out.summary_var_by_batch >> 'summary_var_by_batch'
  VIRAL_VARIANT.out.var_batch_filtered >> 'var_batch_filtered'             //
  VIRAL_VARIANT.out.var_by_smpl_filtered >> 'var_by_smpl_filtered'
  VIRAL_VARIANT.out.var_by_smpl_corrected >> 'var_by_smpl_corrected'
  VIRAL_VARIANT.out.subset_mpileup_cov_by_smpl >> 'subset_mpileup_cov_by_smpl'     //
  VIRAL_VARIANT.out.transfered_gff >> 'transfered_gff'
  VIRAL_VARIANT.out.psa_algn >> 'psa_algn'                        //
  VIRAL_VARIANT.out.psa_genomic_coords >> 'psa_genomic_coords'   //
  VIRAL_VARIANT.out.flagstat >> 'flagstat'
  // bam in option: save_bam ?

}



output {
  directory 'out_dir'
  mode 'rellink'
}

