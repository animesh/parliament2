version 1.0

workflow Parliament2 {
    input {
        File inputBam
        File inputBai
        File refFasta 
        File refIndex
        Boolean filterContigs
        Boolean runBreakdancer
        Boolean runBreakseq
        Boolean runCNVnator
        Boolean runManta
    }

    if (runBreakdancer || runCNVnator || runManta) {
        call P2Prep {
            input:
                inputBam = inputBam,
                filterContigs = filterContigs
        }
        
        if (runCNVnator) {
            Array[String] chromosomes = read_lines(P2Prep.contigs)
            scatter (chromosome in chromosomes) {
                call CNVnatorChrom {
                    input:
                        inputBam = inputBam,
                        inputBai = inputBai,
                        refFasta = refFasta,
                        refIndex = refIndex,
                        contig = chromosome
                }
            }
        }
    }

    # output {
    #     File? breakdancerCTX = GatherBreakdancer.breakdancerCTX
    #     File? breakdancerVCF = GatherBreakdancer.breakdancerVCF

    #     File? breakseqGFF = Breakseq.breakseqGFF
    #     File? breakseq_genotypedGFF = Breakseq.breakseq_genotypedGFF
    #     File? breakseqVCF = Breakseq.breakseqVCF
    #     File? breakseqVCFindex = Breakseq.breakseqVCFindex
    #     File? breakseqBAM = Breakseq.breakseqBAM

    #     File? mantaVCF = Manta.mantaVCF
    #     File? mantaStats = Manta.mantaStats
    #     File? mantaVariants = Manta.mantaVariants
    # }
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

task CNVnatorChrom {
    input {
        File inputBam
        File inputBai
        File refFasta
        File refIndex
        String contig
    }

    command <<<
        cnvnator -root output.root"~{contig}" -chrom "~{contig}" -genome "~{refFasta}" -tree "~{inputBam}"
        cnvnator -root output.root"~{contig}" -chrom "~{contig}" -genome "~{refFasta}" -his 100
        cnvnator -root output.root"~{contig}" -chrom "~{contig}" -genome "~{refFasta}" -stat 100
        cnvnator -root output.root"~{contig}" -chrom "~{contig}" -genome "~{refFasta}" -partition 100
        cnvnator -root output.root"~{contig}" -chrom "~contig" -genome "~{refFasta}" -call 100  > "output.cnvnator_calls-~{contig}"
    >>>

    Int diskGb = ceil(2.0 * size(inputBam, "G"))

    runtime {
        docker : "szarate/cnvnator:v0.4.1"
        disks : "local-disk ${diskGb} SSD"
        cpu : 8
    }

    output {
        File CNVnatorOutput = "output.cnvnator_calls-~{contig}"
    }
}
