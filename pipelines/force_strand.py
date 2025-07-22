def forceStrand(input, output):
    with open(input) as infile, open(output, "w") as outfile:
        for line in infile:
            outfile.write(line.replace('\t?\t', '\t-\t'))

# example use
# forceStrand('/gpfs/commons/home/ncui/project/falsitron_pipeline/na12878_cDNA.stranded.junc', '/gpfs/commons/home/ncui/project/falsitron_pipeline/na12878_cDNA.fixed.stranded.junc')
