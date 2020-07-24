import sys

def main():
    filter_contigs = sys.argv[1]

    for line in sys.stdin:
        if line[:3] != "@SQ": 
            continue
        sequence_name = line.strip().split("SN:")[-1].split("\t")[0]
        length = int(line.strip().split("LN:")[-1].split("\t")[0])

        # Ignore the sponge contig for the hs37d5 human genome
        if filter_contigs and (sequence_name == "hs37d5" or "alt" in sequence_name or "_random" in sequence_name or "_decoy" in sequence_name):
            continue

        if filter_contigs and length < 1000000:
            continue

        if filter_contigs and (length < 1000000 or sequence_name == "hs37d5" or "alt" in sequence_name or "_random" in sequence_name or "_decoy" in sequence_name):
            continue
        else:
            print(sequence_name)

main()