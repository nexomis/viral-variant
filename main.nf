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
  log.info showSchemaHelp("assets/input_schema.json")
  log.info showSchemaHelp("assets/refs_schema.json")
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
    "batch_id": it[3],
    "rank_in_batch": it[4]
  ]
  return [meta, files]
}
  
// include
include {PRIMARY} from './modules/subworkflows/primary/main.nf'
include {VIRAL_VARIANT} from './modules/subworkflows/viral_variant/main.nf'

workflow {

  Channel.fromSamplesheet("input")
  | map {
    return parse_sample_entry(it)
  }
  | set { parsedInputs }

  if (params.skip_primary) {
    inReads = parsedInputs
  } else {
    if ( params.kraken2_db == null ) {
      error "kraken2_db argument required for primary analysis"
    }

    Channel.fromPath(params.kraken2_db, type: "dir", checkIfExists: true)
    | map {[["id": "kraken_db"], it]}
    | collect
    | set {dbPathKraken2}

    taxDir = Channel.fromPath(params.tax_dir, type: 'dir')

    numReads = Channel.value(params.num_reads_sample_qc)
    
    PRIMARY(parsedInputs, dbPathKraken2, taxDir, numReads)
    PRIMARY.out.trimmed
    | set { inReads }
  }

  inputRefs = Channel.fromSamplesheet("refs")

  inRef = inputRefs.map {[[id: it[0]], it[3] == "" ? [file(it[1]), file(it[2])] : [file(it[1])]]}
  inAnnot = inputRefs.filter{it[3] != ""}.map{[[id: it[0]],[file(it[3]), file(it[2])]]}

  VIRAL_VARIANT(inReads, inRef, inAnnot, params.mapper)
  
  publish:

  PRIMARY.out.trimmed >> 'trimmed_and_filtered'
  PRIMARY.out.fastqc_trim_html >> 'fastqc_for_trimmed'
  PRIMARY.out.fastqc_raw_html >> 'fastqc_for_raw'
  PRIMARY.out.multiqc_html >> 'multiqc'
  PRIMARY.out.kraken2_report >> 'classification'
  PRIMARY.out.kraken2_output >> 'classification'
  PRIMARY.out.class_report >> 'classification'

  VIRAL_VARIANT.out.called_snv >> 'sav_call_per_sample/snv'
  VIRAL_VARIANT.out.called_snv_rev >> 'sav_call_per_sample/snv_rev'
  VIRAL_VARIANT.out.called_snv_fwd >> 'sav_call_per_sample/snv_fwd'
  VIRAL_VARIANT.out.called_indel >> 'sav_call_per_sample/indel'
  VIRAL_VARIANT.out.sav_base >> 'sav_call_per_sample/details'

  VIRAL_VARIANT.out.variants >> 'variants_per_batch'
  VIRAL_VARIANT.out.vcf >> 'variants_per_batch'
  VIRAL_VARIANT.out.proteins >> 'variants_per_batch'

  VIRAL_VARIANT.out.transfered_gff >> 'transfered_annot'
  VIRAL_VARIANT.out.psa_align >> 'transfered_annot'
  VIRAL_VARIANT.out.flagstat >> 'mapping/flagstat'
  VIRAL_VARIANT.out.alignedBam >> 'mapping'


}

output {
    'trimmed_and_filtered' {
        enabled = params.trimmed_and_filtered
    }
    'mapping' {
        enabled = params.save_bam_after_indel_realn
    }
}
