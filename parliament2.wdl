workflow parliament2{
    File inputBAM
    File inputBAI

    File fasta
    File fastaFAI

    Int threads

    call breakdancer{
        input:
            inputBAM=inputBAM,
            inputBAI=inputBAI,
            fasta=fasta,
            fastaFAI=fastaFAI,
            threads=threads
    }
}

task breakdancer{
    File inputBAM
    File inputBAI

    File fasta
    File fastaFAI

    Int threads

    String bamBase = basename(inputBAM)
    String baiBase = basename(inputBAI)
    String faBase = basename(fasta)
    String faiBase = basename(fastaFAI)

    command {
        /miniconda/bin/breakdancer-max
        /miniconda/bin/bam2cfg.pl
    }

    runtime {
        docker : "slzarate/parliament2"
        disks : "local-disk 1000 HDD"
        memory : "62G"
        cpu : "32"
    }

    output{
        File lumpyVCF = "lumpy.vcf"
        File cnvnatorVCF = "cnvnator.vcf"
        File lumpyGFF = "lumpy.gff"
        Array[File] dellyVCFs = glob("*.delly.deletion.vcf")
        File breakdancerVCF = "breakdancer.vcf"
    }
}