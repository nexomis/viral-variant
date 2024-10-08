process FILTER_REGROUP_IVAR_VARIANTS {
  container "${params.biocontainers_registry ?: 'quay.io'}/biocontainers/pandas:2.2.1"   // 117 MiB !!

  input:
  tuple val(meta), path(ivar_tsv_files, arity: '1..*', stageAs: 'input/')  // 'stageAs' and templates ?
  tuple val(meta2), path(pos_file, arity: 1, stageAs: 'input/')  // 'stageAs' and templates ?

  output:
  tuple val(meta), path("${meta.id}_summary_all_iSNVs.tsv", type: 'file') , optional:false ,  emit: batch_summary_all_iSNVs
  tuple val(meta), path("${meta.id}_batchFiltered", type: 'dir') , optional:false ,  emit: smpl_summary_all_iSNVs

  script:
  //id 'pos_file' is empty, this occure bug, this is improbable but is better to manage it.
  template "filter_regroup_ivar_variants.py"

  stub:
  """
  #!/usr/bin/bash

  touch ${meta.id}_summary_all_iSNVs.tsv
  mkdir ${meta.id}/
  """
}
