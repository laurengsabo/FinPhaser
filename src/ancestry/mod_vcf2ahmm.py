import argparse
import gzip
import csv
from collections import defaultdict

### get the args
parser = argparse.ArgumentParser()
parser.add_argument("-v", type=str, help="input vcf file", required=True )
parser.add_argument("-s", type=str, help="sample to population file", required=True )
parser.add_argument("-o_txt", type=str, help="output text file name", default="ancestry_input.txt")
parser.add_argument("-o_ploidy", type=str, help="ploidy file for ahmm input", default = "ahmm.ploidy")
parser.add_argument("-g", type=int, help="boolean (0 = use sample reads| 1 = use genotypes)", default = 0 )
parser.add_argument("-r", type=float, help="recombination rate (Morgans/bp)", default = 1e-8 )
parser.add_argument("-m", type=int, help="min distance between SNPs (bp)", default = 1000 )
parser.add_argument("--min_total", type=int, help="min number of reference alleles", default = 1 )
parser.add_argument("--min_diff", type=float, help="min allele freq difference", default = 0.1 )
args = parser.parse_args()

### Load mappings
sample2pop = defaultdict(list)
with open(args.s) as file:
    for line in file:
        parts = line.split()
        if len(parts) >= 2:
            sample2pop[parts[0]].append(parts[1])

# Clear the ploidy file if it exists so we don't append to old data
with open(args.o_ploidy, "w") as pf:
    pass

sample_names = []
last_position = -10000000
current_chrom = "NA"
ploidy_print = 0

# Open the main output file
with open(args.v) as tsv, open(args.o_txt, "w") as out_f:
    for line in csv.reader(tsv, delimiter="\t"):
        if line[0].startswith('#CHROM'):
            sample_names = line
            continue
        if line[0].startswith('##'):
            continue
        
        if len(line[3]) != 1 or len(line[4]) != 1:
            continue

        pos = int(line[1])
        if current_chrom != line[0]:
            current_chrom = line[0]
            last_position = -1000000
        elif pos - last_position < args.m:
            continue

        alts = defaultdict(int)
        totals = defaultdict(int)
        
        for s in range(9, len(line)):
            name = sample_names[s]
            if name in sample2pop:
                pops = sample2pop[name]
                if "admixed" in pops:
                    continue
                
                gt_field = line[s].split(':')[0]
                if '|' not in gt_field and '/' in gt_field:
                    continue
                    
                haplotypes = gt_field.split('|')
                for i, allele in enumerate(haplotypes):
                    if i < len(pops):
                        p_id = pops[i]
                        if allele != '.':
                            totals[p_id] += 1
                            if allele == '1':
                                alts[p_id] += 1

        sorted_pops = sorted([p for p in totals.keys()])
        
        # Ensure all 4 haplotypes (0, 1, 2, 3) are represented
        if len(sorted_pops) < 4:
            continue

        freqs = [alts[p]/totals[p] if totals[p] > 0 else 0 for p in sorted_pops]
        diff_found = False
        for i in range(len(freqs)):
            for j in range(i + 1, len(freqs)):
                if abs(freqs[i] - freqs[j]) > args.min_diff:
                    diff_found = True
                    break
        
        if not diff_found:
            continue

        # Write to the .txt file instead of printing to console
        out_f.write(f"{line[0]}\t{line[1]}\t")
        
        for p in sorted_pops:
            ref_c = totals[p] - alts[p]
            out_f.write(f"{ref_c}\t{alts[p]}\t")

        rec = (pos - last_position) * args.r if last_position > 0 else 0
        out_f.write(f"{rec:.10f}")

        for s in range(9, len(line)):
            name = sample_names[s]
            if name in sample2pop and "admixed" in sample2pop[name]:
                if args.g == 1:
                    gt = line[s].split(':')[0].replace('|', '/')
                    gt_alleles = gt.split('/')
                    out_f.write(f"\t{gt_alleles.count('0')}\t{gt_alleles.count('1')}")
                else:
                    fields = line[s].split(':')
                    ad = fields[1].split(',') if len(fields) > 1 else ["0", "0"]
                    out_f.write(f"\t{ad[0]}\t{ad[1]}")
                
                if ploidy_print == 0:
                    with open(args.o_ploidy, "a") as pf:
                        pf.write(f"{name}\t2\n")

        out_f.write("\n")
        last_position = pos
        ploidy_print = 1

print(f"Done! Created {args.o_txt} and {args.o_ploidy}")