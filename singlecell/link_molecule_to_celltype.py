from csv import DictReader, DictWriter
import gzip
import os, sys

def main(bc_info_csv, input_csv, output_csv):
    bc_info = {}  # BC --> cell type
    reader = DictReader(open(bc_info_csv), delimiter=',')
    if 'BCrev' not in reader.fieldnames or 'Celltype' not in reader.fieldnames:
        print(f"ERROR: Expected 'BCrev' and 'Celltype' columns in {bc_info_csv}. Abort!")
        sys.exit(-1)
    for r in reader:
        bc_info[r['BCrev']] = r['Celltype']

    reader = DictReader(gzip.open(input_csv, 'rt'), delimiter='\t')
    if 'BCrev' not in reader.fieldnames or 'id' not in reader.fieldnames:
        print(f"ERROR: Expected 'id' and 'BCrev' columns in {input_csv}. Abort!")
        sys.exit(-1)

    f = open(output_csv, 'w')
    writer = DictWriter(f, fieldnames=reader.fieldnames+['Celltype'], delimiter=',')
    writer.writeheader()
    for r in reader:
        if r['BCrev'] in bc_info:
            r['Celltype'] = bc_info[r['BCrev']]
            writer.writerow(r)
    f.close()
    print(f"Output written to: {output_csv}")

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("bc_csv", help="Comma-delimited barcode to cell type information, must have 'BCrev' and 'Celltype' columns")
    parser.add_argument("annotated_csv", help="Tab-delimited annotated CSV file, gzipped, must have 'id' and 'BCrev' columns")
    parser.add_argument("output_csv", help="Output CSV filename")

    args = parser.parse_args()
    if not os.path.exists(args.bc_csv):
        print("Input file {0} does not exist! Abort!".format(args.bc_csv))
        sys.exit(-1)
    if not os.path.exists(args.annotated_csv):
        print("Input file {0} does not exist! Abort!".format(args.annotated_csv))
        sys.exit(-1)

    main(args.bc_csv, args.annotated_csv, args.output_csv)