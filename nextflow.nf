$HOSTNAME = ""
params.outdir = 'results'  



if (!params.inputs){params.inputs = ""} 
if (!params.vector_fa){params.vector_fa = ""} 
if (!params.packaging_fa){params.packaging_fa = ""} 
if (!params.host_fa){params.host_fa = ""} 
if (!params.mate){params.mate = ""} 
if (!params.flipflop_fa){params.flipflop_fa = ""} 
if (!params.vector_bed){params.vector_bed = ""} 
if (!params.bed_fn){params.bed_fn = ""} 
if (!params.vcf_fn){params.vcf_fn = ""} 
// Stage empty file to be used as an optional input where required
ch_empty_file_1 = file("$baseDir/.emptyfiles/NO_FILE_1", hidden:true)
ch_empty_file_2 = file("$baseDir/.emptyfiles/NO_FILE_2", hidden:true)

if (params.inputs){
Channel
	.fromFilePairs( params.inputs,checkExists:true , size: params.mate == "single" ? 1 : params.mate == "pair" ? 2 : params.mate == "triple" ? 3 : params.mate == "quadruple" ? 4 : -1 ) 
	.set{g_4_0_g_1}
  } else {  
	g_4_0_g_1 = Channel.empty()
 }

g_5_1_g_1 = file(params.vector_fa, type: 'any')
g_5_1_g_13 = file(params.vector_fa, type: 'any')
g_6_2_g_1 = params.packaging_fa && file(params.packaging_fa, type: 'any').exists() ? file(params.packaging_fa, type: 'any') : ch_empty_file_1
g_7_3_g_1 = params.host_fa && file(params.host_fa, type: 'any').exists() ? file(params.host_fa, type: 'any') : ch_empty_file_2
Channel.value(params.mate).set{g_8_4_g_1}
g_9_1_g_2 = params.flipflop_fa && file(params.flipflop_fa, type: 'any').exists() ? file(params.flipflop_fa, type: 'any') : ch_empty_file_1
g_10_2_g_2 = params.vector_bed && file(params.vector_bed, type: 'any').exists() ? file(params.vector_bed, type: 'any') : ch_empty_file_2
g_15_2_g_13 = params.bed_fn && file(params.bed_fn, type: 'any').exists() ? file(params.bed_fn, type: 'any') : ch_empty_file_1
g_16_3_g_13 = params.vcf_fn && file(params.vcf_fn, type: 'any').exists() ? file(params.vcf_fn, type: 'any') : ch_empty_file_2


process map_reads {

input:
 tuple val(name), file(inputs)
 path vector_fa
 path packaging_fa
 path host_fa
 val mate

output:
 tuple val(name), file("${name}.reference_names.tsv"), file("${name}.sort_by_name.bam")  ,emit:g_1_bamFile00_g_2 
 tuple val(name), file("${name}.sort_by_pos.bam"), file("${name}.sort_by_pos.bam.bai")  ,emit:g_1_bam_bai10_g_13 


container 'quay.io/ummsbiocore/laava:1.0.0'

script:

file1 =  inputs[0].toString() 
file2 = ""
if (mate == "pair") {file2 =  " $inputs[1]" }
gunzipCmd = ""
if (file1.endsWith(".bam.gz")) {
    gunzipCmd = "gunzip ${file1}"
    file1 = file1.replaceFirst(/\.gz$/, "")
}
packaging_fa_path = packaging_fa.name.startsWith('NO_FILE') ? "" : "${packaging_fa}"
host_fa_path = host_fa.name.startsWith('NO_FILE') ? "" : "${host_fa}"

"""
$gunzipCmd 
map_reads.sh ${name} '${file1}${file2}' '${vector_fa}' '${packaging_fa_path}' '${host_fa_path}' '${params.repcap_name}' '${params.helper_name}' '${params.lambda_name}'
"""
}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 10
    $MEMORY = 30
}
//* platform
//* platform
//* autofill

process make_report {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${sample_id}_AAV_report.html$/) "html_report/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${sample_id}_AAV_report.pdf$/) "pdf_report/$filename"}
input:
 tuple val(sample_id),  file(reference_names), file(mapped_reads)
 path flipflop_fa
 path vector_annotation

output:
 val sample_id  ,emit:g_2_sample_id00 
 path "${sample_id}.metadata.tsv"  ,emit:g_2_outFileTSV11 
 path "${sample_id}.alignments.tsv.gz"  ,emit:g_2_outFileTSV22 
 path "${sample_id}.per_read.tsv.gz"  ,emit:g_2_outFileTSV33 
 path "${sample_id}.nonmatch.tsv.gz"  ,emit:g_2_outFileTSV44 
 path "${sample_id}.agg_ref_type.tsv"  ,emit:g_2_outFileTSV55 
 path "${sample_id}.agg_subtype.tsv"  ,emit:g_2_outFileTSV66 
 path "${sample_id}.flipflop.tsv.gz" ,optional:true  ,emit:g_2_outFileTSV77 
 path "${sample_id}.agg_flipflop.tsv" ,optional:true  ,emit:g_2_outFileTSV88 
 path "${sample_id}.tagged.bam"  ,emit:g_2_bamFile99 
 path "${sample_id}.*.tagged.sorted.bam"  ,emit:g_2_bamFile1010 
 path "${sample_id}.*.tagged.sorted.bam.bai"  ,emit:g_2_bam_bai1111 
 path "${sample_id}.flipflop-*.bam" ,optional:true  ,emit:g_2_bamFile1212 
 path "${sample_id}_AAV_report.html"  ,emit:g_2_outputFileHTML1313 
 path "${sample_id}_AAV_report.pdf"  ,emit:g_2_outputFilePdf1414 

container 'quay.io/ummsbiocore/laava:1.0.0'

script:
ff_fa_path = flipflop_fa.name.startsWith('NO_FILE') ? "" : "${flipflop_fa}"
"""

name=\$(basename "${vector_annotation}")

if [[ "\$name" == NO_FILE* ]]; then
	
    IFS='-' read -r start1 end1 <<< "${params.ITR_pos_1}"
    IFS='-' read -r start2 end2 <<< "${params.ITR_pos_2}"
    IFS='-' read -r start3 end3 <<< "${params.mITR_pos}"

    # write the two lines
    {
      if [[ -n \$start1 && -n \$end1 ]]; then
    	printf "pAV_CMV_GFP\\t\$start1\\t\$end1\\tITR\\n"
      fi
      if [[ ${params.vector_type} == "ss" && -n \$start2 && -n \$end2 ]]; then
    	printf "pAV_CMV_GFP\\t\$start2\\t\$end2\\tITR\\n"
      fi
      if [[ ${params.vector_type} == "sc" && -n \$start3 && -n \$end3 ]]; then
    	printf "pAV_CMV_GFP\\t\$start3\\t\$end3\\tmITR\\n"
      fi
    } > 'custom.bed'
    final_bed='custom.bed'
    itr_lbl_1="ITR"
    itr_lbl_2="ITR"
    mitr_lbl="mITR"
else
	final_bed=${vector_annotation}
	itr_lbl_1=${params.itr_label_1}
	itr_lbl_2=${params.itr_label_2}
	mitr_lbl=${params.mitr_label}
fi

echo \$name
echo \$final_bed
head \$final_bed

make_report.sh \
        '${sample_id}' \
        '${sample_id}' \
        '${params.manifest_version}' \
        '${reference_names}' \
        '${mapped_reads}' \
        "\$final_bed" \
        "\$itr_lbl_1" \
        "\$itr_lbl_2" \
        "\$mitr_lbl" \
        '${params.vector_type}' \
        '${params.target_gap_threshold}' \
        '${params.max_allowed_outside_vector}' \
        '${params.max_allowed_missing_flanking}' \
        '${params.min_supp_joint_coverage}' \
        '${params.flipflop_name}' \
        '${ff_fa_path}' \
        '.'

"""
}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 30
    $MEMORY = 20
}
//* platform
//* platform
//* autofill

process Clair3 {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${sample_name}.vcf.gz$/) "variants/$filename"}
input:
 tuple val(sample_name), file(input_bam), file(input_bai)
 path vector_fa
 path bed_fn
 path vcf_fn

output:
 path "${sample_name}.vcf.gz"  ,emit:g_13_outputFileVCF00 

container "hkubal/clair3:v1.0.11"

when:
(params.run_Clair3 && (params.run_Clair3 == "yes")) || !params.run_Clair3

script:

platform = params.Clair3.platform
model_name = params.Clair3.model_name
threads = task.cpus

opt_params = params.Clair3.opt_params
bed_fn_path = bed_fn.name.startsWith('NO_FILE') ? "" : "--bed_fn=${bed_fn}"
vcf_fn_path = vcf_fn.name.startsWith('NO_FILE') ? "" : "--vcf_fn=${vcf_fn}"

"""

mkdir -p ${sample_name}

samtools faidx ${vector_fa}
	
/opt/bin/run_clair3.sh \
  --bam_fn='${input_bam}' \
  --ref_fn=${vector_fa} \
  --threads='${threads}' \
  --platform=${platform} \
  --model_path='/opt/models/${model_name}' \
  --output='${sample_name}' \
  --include_all_ctgs \
  ${opt_params} ${bed_fn_path} ${vcf_fn_path}

mv ${sample_name}/merge_output.vcf.gz ${sample_name}.vcf.gz
"""

}


workflow {



map_reads(g_4_0_g_1,g_5_1_g_1,g_6_2_g_1,g_7_3_g_1,g_8_4_g_1)
g_1_bamFile00_g_2 = map_reads.out.g_1_bamFile00_g_2
g_1_bam_bai10_g_13 = map_reads.out.g_1_bam_bai10_g_13



make_report(g_1_bamFile00_g_2,g_9_1_g_2,g_10_2_g_2)
g_2_sample_id00 = make_report.out.g_2_sample_id00
g_2_outFileTSV11 = make_report.out.g_2_outFileTSV11
g_2_outFileTSV22 = make_report.out.g_2_outFileTSV22
g_2_outFileTSV33 = make_report.out.g_2_outFileTSV33
g_2_outFileTSV44 = make_report.out.g_2_outFileTSV44
g_2_outFileTSV55 = make_report.out.g_2_outFileTSV55
g_2_outFileTSV66 = make_report.out.g_2_outFileTSV66
g_2_outFileTSV77 = make_report.out.g_2_outFileTSV77
g_2_outFileTSV88 = make_report.out.g_2_outFileTSV88
g_2_bamFile99 = make_report.out.g_2_bamFile99
g_2_bamFile1010 = make_report.out.g_2_bamFile1010
g_2_bam_bai1111 = make_report.out.g_2_bam_bai1111
g_2_bamFile1212 = make_report.out.g_2_bamFile1212
g_2_outputFileHTML1313 = make_report.out.g_2_outputFileHTML1313
g_2_outputFilePdf1414 = make_report.out.g_2_outputFilePdf1414



Clair3(g_1_bam_bai10_g_13,g_5_1_g_13,g_15_2_g_13,g_16_3_g_13)
g_13_outputFileVCF00 = Clair3.out.g_13_outputFileVCF00


}

workflow.onComplete {
println "##Pipeline execution summary##"
println "---------------------------"
println "##Completed at: $workflow.complete"
println "##Duration: ${workflow.duration}"
println "##Success: ${workflow.success ? 'OK' : 'failed' }"
println "##Exit status: ${workflow.exitStatus}"
}
