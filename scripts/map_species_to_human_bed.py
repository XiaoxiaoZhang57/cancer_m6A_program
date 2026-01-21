import csv


def parse_paf(paf_file):
    mappings = {}
    with open(paf_file, 'r') as f:
        for line in f:
            fields = line.strip().split()
            qname = fields[0]  # Query  chromosome
            qstart = int(fields[2])  # Query start position
            qend = int(fields[3])  # Query end position
            tname = fields[5]  # Target (human) chromosome
            tstart = int(fields[7])  # Target start position
            tend = int(fields[8])  # Target end position
            mappings[qname] = mappings.get(qname, []) + [(qstart, qend, tname, tstart, tend)]
    return mappings


def convert_positions(bed_file, mappings, output_file):
    with open(bed_file, 'r') as bed, open(output_file, 'w', newline='') as out:
        bed_reader = csv.reader(bed, delimiter='\t')
        out_writer = csv.writer(out, delimiter='\t')
        for row in bed_reader:
            achr, astart, aend = row[0], int(row[1]), int(row[2])
            if achr in mappings:
                for qstart, qend, hchr, hstart, hend in mappings[achr]:
                    if qstart <= astart <= qend or qstart <= aend <= qend:
                    
                        offset_start = astart - qstart
                        offset_end = aend - qstart
                        human_start = hstart + offset_start
                        human_end = hstart + offset_end
                        out_writer.writerow([achr, astart, aend, hchr, human_start, human_end])
                        break


paf_file = 'A_to_human.paf'
bed_file = 'A.bed'
output_file = 'A_to_human.bed'
mappings = parse_paf(paf_file)
convert_positions(bed_file, mappings, output_file)


