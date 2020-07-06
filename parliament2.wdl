version 1.0

workflow Parliament2{
    input {
        File inputBam
        File inputBai
        File refFasta 
        File refIndex
        Boolean runBreakseq
    }

    if (runBreakseq) {
        call Breakseq { 
            input: 
                inputBam = inputBam,
                inputBai = inputBai,
                refFasta = refFasta,
                refIndex = refIndex
        }
    }

}

task Breakseq {
    input {
        File inputBam
        File inputBai
        File refFasta
        File refIndex
    }

    String bamBase='~{basename(inputBam,".bam")}'
    String refBase = basename(refFasta)
    String breakpointLibrary = if (refBase == "*hg19*") then "/breakseq2_bplib_20150129.hg19/breakseq2_bplib_20150129.hg19.gff" else (if refBase == "*hg38*" then "/bplib.hg38.gff" else "/breakseq2_bplib_20150129.hs37d5/breakseq2_bplib_20150129.gff")

    command <<<
        mkdir -p "breakseq2"
        gunzip "~{refFasta}"

        refName="~{refFasta}"
        refName="${refName%.gz}"

        /miniconda/bin/run_breakseq2.py \
            --reference "${refName}" \
            --bams "~{inputBam}" \
            --work breakseq2 \
            --bwa /miniconda/bin/bwa \
            --samtools /miniconda/bin/samtools \
            --bplib_gff "~{breakpointLibrary}" \
            --nthreads "$(nproc)" \
            --sample "~{bamBase}"

            mv breakseq2/breakseq.vcf.gz "~{bamBase}.vcf.gz"
            mv breakseq2/breakseq.vcf.gz.tbi "~{bamBase}.vcf.gz.tbi"
            mv breakseq2/breakseq.gff "~{bamBase}.gff"
            mv breakseq2/breakseq_genotyped.gff "~{bamBase}_genotyped.gff"
            mv breakseq2/final.bam "~{bamBase}.breakseq.bam"
    >>>

    runtime {
        docker : "szarate/breakseq2:v2.2"
    }

    output {
        File breakseqGFF = "${bamBase}.gff"
        File breakseq_genotypedGFF = "${bamBase}_genotyped.gff"
        File breakseqVCF = "${bamBase}.vcf.gz"
        File breakseqVCFindex = "${bamBase}.vcf.gz.tbi"
        File breakseqBAM = "${bamBase}.breakseq.bam"
    }
}