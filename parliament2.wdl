version 1.0

workflow Parliament2 {
    input {
        File inputBam
        File inputBai
        File refFasta 
        File refIndex
        Boolean filterContigs
        Boolean runBreakseq
        Boolean runManta
    }

    if (runManta) {
        call P2Prep {
            input:
                inputBam = inputBam,
                filterContigs = filterContigs
        }

        call Manta {
            input:
                inputBam = inputBam,
                inputBai = inputBai,
                refFasta = refFasta,
                refIndex = refIndex,
                contigs = P2Prep.contigs
        }
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

task P2Prep {
    input {
        File inputBam
        Boolean filterContigs
    }

    command <<<
        samtools view -H "~{inputBam}" | python /get_contigs.py "~{filterContigs}" > contigs
    >>>

    Int diskGb = ceil(2.0 * size(inputBam, "G"))
    
    runtime {
        docker : "szarate/p2_prep:v0.0.1"
        disks : "local-disk ${diskGb} SSD"
    }

    output {
        File contigs = "contigs"
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
        File breakseqGFF = "${bamBase}.breakseq.gff"
        File breakseq_genotypedGFF = "${bamBase}_genotyped.breakseq.gff"
        File breakseqVCF = "${bamBase}.breakseq.vcf.gz"
        File breakseqVCFindex = "${bamBase}.breakseq.vcf.gz.tbi"
        File breakseqBAM = "${bamBase}.breakseq.bam"
    }
}

task Manta {
    input {
        File inputBam
        File inputBai
        File refFasta
        File refIndex
        File contigs
    }

    String bamBase='~{basename(inputBam,".bam")}'

    command <<<
        mkdir -p "manta"
        gunzip "~{refFasta}"

        refName="~{refFasta}"
        refName="${refName%.gz}"
        region_string=""

        while read line; do
            region_string="$region_string --region=$line"
        done < "~{contigs}"

        python /usr/local/bin/configManta.py --referenceFasta "${refName}" --normalBam "~{inputBam}" --runDir manta $region_string

        python manta/runWorkflow.py -m local -j "$(nproc)"

        tar -czf stats.tar.gz -C manta/results/stats/ .
        tar -czf variants.tar.gz -C manta/results/variants/ .

        mv manta/results/variants/diploidSV.vcf.gz "~{bamBase}.manta.vcf.gz"
        mv stats.tar.gz "~{bamBase}_stats.manta.vcf.gz"
        mv variants.tar.gz "~{bamBase}_variants.manta.vcf.gz"
    >>>

    runtime {
        docker : "szarate/manta:v1.6.0"
    }

    output {
        File mantaVCF = "~{bamBase}.manta.vcf.gz"
        File mantaStats = "~{bamBase}_stats.manta.vcf.gz"
        File mantaVariants = "~{bamBase}_variants.manta.vcf.gz"
    }
}