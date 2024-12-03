from Bio import SeqIO

def extract_locus_tag_protein_id(genbank_file, output_file):
    # Open the output CSV file
    with open(output_file, "w") as out_file:
        # Write the header
        out_file.write("locus_tag,protein_id\n")
        
        # Parse the GenBank file
        for record in SeqIO.parse(genbank_file, "genbank"):
            # Iterate through features to find CDS features
            for feature in record.features:
                if feature.type == "CDS":
                    # Extract the locus_tag
                    locus_tag = feature.qualifiers.get("locus_tag", [None])[0]
                    
                    # Extract the protein_id
                    protein_id = feature.qualifiers.get("protein_id", [None])[0]
                    
                    # Only write to file if both locus_tag and protein_id are present
                    if locus_tag and protein_id:
                        out_file.write(f"{locus_tag},{protein_id}\n")

# Usage example
extract_locus_tag_protein_id("NC_007880.gb", "NC_007880_output_locus_refseq.csv")
