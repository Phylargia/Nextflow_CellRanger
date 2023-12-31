#!/usr/bin/env nextflow

// Help Page Documentation
if( params.help ) {

log.info """
=====================================
Help Documentation : "Exome Sequencing Nextflow Pipeline"

Developed by Christopher O'Connor

Current Version: Pre-Alpha 
=====================================
Usage:
          Exome Sequencing of fastq reads using the cellranger module 
	  Example: nextflow run main_new.nf --genome mouse --input 'input/input.csv'

Arguments:
          * --input  : Option to summon an array matrix defining sample(s)/batch provided with pathway  
          * --genome : Option to define define transcriptome reference (Human/GRCh38, Mouse/mm10) 
          * --help   : Option to call on help documentation page to view all available arguments 
=====================================
Configuration:
* In nextflow.config, please adjust the parameters for the following: *
	  
      	  ** Priority List **
      	  * <userID>          : Define scw user ID
          * <userProjectCode> : Define scw project ID  
          * <cellRanger>      : Define pathway to cellranger module 
      	  * <projectID>       : Define project name 
      	  * <dataID>          : Define output data name
=====================================
Processes:
          1) Read input file 
      	  2) FastQC            
      	  3) MultiQC
      	  4) Cellranger Count
      	  5) Cellranger Aggregate
=====================================
Outputs:
      	  * FastQC reports 
      		└── Quality control check on raw sequence data	  
      
      	  * MultiQC report (HTML)
      		└── Aggregated report of FastQC reports
	   
      	  * Cellranger Count / Cellranger Aggregate
                └── Aggregation metrics summary HTML       : ../outs/web_summary.html
                └── Copy of the input aggregation CSV      : ../outs/aggregation.csv
                └── Aggregation metrics summary JSON 	   : ../outs/summary.json
                └── Secondary analysis output CSV          : ../outs/analysis 
                └── Filtered feature-barcode matrices MEX  : ../outs/filtered_feature_bc_matrix 
            	└── Filtered feature-barcode matrices HDF5 : ../outs/filtered_feature_bc_matrix.h5
            	└── Loupe Browser file                     : ../outs/cloupe.cloupe
=====================================
Notes:
      	* Pipeline remains in early development (Please contact kit.ourea@gmail.com for assistance)
      	* Pipeline has been tested and validated to perform with fastq read files approx. 50GB in size	
        * Pipeline has been designed to run using a 'slurm' operator 
        * Error strategy has been designed to repeat once on failure
      	* Pipeline has been designed to be followed up with scRNAseq analysis using Seurat or Scanpy

Roadmap: 
        * Changes to '--genome' to allow custom references 
      	* Improve error handling
      	* Improve directory handling
      	* Potentially, consider follow-up R/Py scripts to perform 'QC' and 'Pre-processing' scRNAseq analysis   
=====================================
"""
        exit 0
} 


// Define Channels
if (params.input)  { ch_metadata = file(params.input, checkIfExists: true) } else { exit 1, "Please provide input file with sample metadata with the '--input' option." }
    Channel.from( ch_metadata )
            .splitCsv(header: true, sep:'\t')
            .map { col -> tuple("${col.Sample}", file("${col.R1}", checkifExists: true),file("${col.R2}", checkifExists: true)) }
            .dump()
            .into{ ch_read_fastqc; ch_read_count }


/* 
   --------------------------------------------------------------------
   ------------------------- Start of script --------------------------
   --------------------------------------------------------------------
*/ 

// Run FastQC on input fastq files specified in '--input' argument 

process fastQC {

    errorStrategy { task.exitStatus != 0 ? 'retry' : 'terminate' }
    maxRetries params.retries
    maxErrors -1
    cache 'lenient'

    tag "$Sample"
    publishDir path:{params.outputDir},mode: params.publish_dir_mode

    module params.fastqcVersion

input:
    tuple val(Sample), file(R1), file(R2) from ch_read_fastqc

output:
    file ("${Sample}_fastqc/*.zip") into fastqc_files

script:
    """
    sleep ${params.sleepTimeStart}
    
    mkdir ${Sample}_fastqc 
    fastqc -o ${Sample}_fastqc -t 2 ${params.dataDir}/fastq/${Sample}/${R1} ${params.dataDir}/fastq/${Sample}/${R2} 

    sleep ${params.sleepTimeEnd}
    """
}



// Run MultiQC on FastQC files 

process multiQC {

    errorStrategy { task.exitStatus != 0 ? 'retry' : 'terminate' }
    maxRetries params.retries
    maxErrors -1
    cache 'lenient'

    tag "Multi"
    publishDir path:{params.outputDir},mode: params.publish_dir_mode

    module params.multiqcVersion
    
input:
    file('*') from fastqc_files.collect()

output:
    file ("multiqc_report.html")

script:
    """
    sleep ${params.sleepTimeStart}

    multiqc .

    sleep ${params.sleepTimeEnd}
    """
}



// Run CellRanger Count

process cellranger_count {

    errorStrategy { task.exitStatus != 0 ? 'retry' : 'terminate' }
    maxRetries params.retries
    maxErrors -1
    cache 'lenient'
    
    tag "$Sample"
    publishDir path:{params.outputDir},mode: params.publish_dir_mode

input: 
    tuple val(Sample), file(R1), file(R2) from ch_read_count
    path transcriptome from params.transcriptome

output:
    file ("${Sample}") into count_files 
    //file ("${Sample}/outs/*.molecule_info.h5") into samplesheet_ch

def ref = ''

if (params.genome) {
  // Check for 'GRCh38' or 'Human'
  if (params.genome == 'GRCh38' || params.genome == 'Human' || params.genome == 'human') {
    ref = 'refdata-gex-GRCh38-2020-A'
  } else if (params.genome == 'mm10' || params.genome == 'Mouse' || params.genome == 'mouse') {
    ref = 'refdata-gex-mm10-2020-A'
  } else {
    print('Invalid genome specified. Please provide a valid genome.')
    return
  }
} else {
  print('Genome is missing. Please provide a valid genome.')
  return
}

script:
    """
    sleep ${params.sleepTimeStart}

    export ${params.cellRanger}/bin
	
    cellranger count \
    --id='${Sample}' \
    --fastqs=${params.dataDir}/fastq/${Sample} \
    --sample='${Sample}' \
    --transcriptome=${transcriptome}/$ref \
    --localcores=4 

    sleep ${params.sleepTimeEnd}
    """
}



// Create an reference array with cellranger count files and pathway 

process aggregation_matrix {

    errorStrategy { task.exitStatus != 0 ? 'retry' : 'terminate' }
    maxRetries params.retries
    maxErrors -1
    cache 'lenient'

    tag "Aggregation"  
    publishDir path:{params.outputDir},mode: params.publish_dir_mode

input:
    file Sample from count_files

output:
    file ("aggr_input_${params.dataID}.csv") into ch_aggr

script:
    """
    sleep ${params.sleepTimeStart}

  if [[ ! -e samplesheet.csv ]]; then
      echo "sample_id,molecule_h5" > aggr_input_${params.dataID}.csv
  fi
     echo "${Sample},${params.projectDir}/${params.outputDir}/${Sample}/outs/molecule_info.h5" >> aggr_input_${params.dataID}.csv

    sleep ${params.sleepTimeEnd}
    """
}



// Run CellRanger Aggregate on CellRanger Count files

process cellranger_aggregate {

    errorStrategy { task.exitStatus != 0 ? 'retry' : 'terminate' }
    maxRetries params.retries
    maxErrors -1
    cache 'lenient'

    tag "Aggregate"
    publishDir path:{params.outputDir},mode: params.publish_dir_mode

input:
    file aggr_input from ch_aggr

output:
    file ("aggregated_${params.dataID}")

script:
    """
    sleep ${params.sleepTimeStart}

    export ${params.cellRanger}/bin

    cellranger aggr \
    --id='aggregated_${params.dataID}' \
    --csv=${aggr_input}

    sleep ${params.sleepTimeEnd}
    """
}

/*
  --------------------------------------------------------------------
  ------------------------- End of script ----------------------------
  --------------------------------------------------------------------
*/

workflow.onComplete {

    println ( workflow.success ? """
        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """ : """
        Failed: ${workflow.errorReport}
        exit status : ${workflow.exitStatus}
        """
    )
}
