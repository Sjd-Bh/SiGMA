import pysam
from matplotlib_venn import venn2
import matplotlib.pyplot as plt

# Function to extract mutations from a VCF file
def get_mutation_positions(vcf_file):
    mutations = set()
    vcf = pysam.VariantFile(vcf_file)
    for record in vcf:
        # Use (chromosome, position) as unique identifier of mutations
        mutations.add((record.chrom, record.pos))
    return mutations

# Load the two VCF files
vcf1_file = "../CellGenome/200kb/snvs_sorted1.vcf.gz"
vcf2_file = "../cellAmpOutput/200/PTA/sim1/subset_sort_rg.vcf"

# Extract mutations from both files
mutations_vcf1 = get_mutation_positions(vcf1_file)
mutations_vcf2 = get_mutation_positions(vcf2_file)

# Print the number of mutations in each VCF
print(f"Number of mutations in {vcf1_file}: {len(mutations_vcf1)}")
print(f"Number of mutations in {vcf2_file}: {len(mutations_vcf2)}")

# Create a Venn diagram to compare the two sets of mutations
venn = venn2([mutations_vcf1, mutations_vcf2], set_labels=('simulated SNVs', 'prosolo detection'))
# Set the title of the plot
plt.title('detected SNVs by prosolo')
# Save the plot as a PNG file
plt.savefig('venn_diagram1.png')

# Display the plot
plt.show()
