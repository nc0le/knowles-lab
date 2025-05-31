import subprocess
import pandas as pd

# Absolute paths to your files
vcf = "/Users/nicolecui/Nicole/CS/KnowlesLab/regtools/tests/integration-test/data/vcf/test1.vcf"
bam = "/Users/nicolecui/Nicole/CS/KnowlesLab/regtools/tests/integration-test/data/bam/test_hcc1395.2.bam"
fasta = "/Users/nicolecui/Nicole/CS/KnowlesLab/regtools/tests/integration-test/data/fa/test_chr22.fa"
gtf = "/Users/nicolecui/Nicole/CS/KnowlesLab/regtools/tests/integration-test/data/gtf/test_ensemble_chr22.2.gtf"

# Run regtools command
command = [
    "regtools", "cis-splice-effects", "identify",
    "-s", "RF",
    vcf, bam, fasta, gtf
]

result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)

# Check if regtools executed properly
if result.returncode != 0:
    print("Regtools failed to execute.")
    print(result.stdout)
    exit(1)

lines = result.stdout.splitlines()

# Find header line
for i, line in enumerate(lines):
    if line.startswith("chrom"):
        header_line = line.strip()
        data_lines = lines[i+1:]
        break
else:
    raise ValueError("Could not find table header starting with 'chrom'.")

columns = header_line.split()
data = []

for line in data_lines:
    if line.strip() == "":
        continue
    tokens = line.strip().split()
    # Only parse rows starting with chromosome numbers
    if not tokens[0].isdigit():
        continue
    data.append(tokens)

# Build pandas dataframe
df = pd.DataFrame(data, columns=columns)

# Done!
print("\nParsed regtools output:")
print(df)
