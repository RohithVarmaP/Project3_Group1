# Open the input and output files
with open('C:/Users/Oliver/Desktop/UoN/Group_Projects/Project3/All_except_big_vcf/All_except_big_vcf/r$
    # Loop through each line in the input file
    for line in infile:
        # Skip comment lines and blank lines
        if line.startswith('#') or line.strip() == '':
            continue

        # Split the line into fields
        fields = line.strip().split('\t')

        # Extract the gene name, chromosome/scaffold name, start position, stop position, and direction
        gene_name = fields[8].split(';')[0].split('=')[1]
        chrom_name = fields[0]
        start_pos = fields[3]
        stop_pos = fields[4]
        direction = fields[6]

        # Write the information to the output file
        outfile.write(f"{gene_name}\t{chrom_name}\t{start_pos}\t{stop_pos}\t{direction}\n")
