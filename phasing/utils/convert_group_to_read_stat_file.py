"""
Converting collapsed group.txt file into a .read_stat file that could be used to run IsoPhase later
"""

import os, sys

def main(group_filename, output_filename):
    f = open(output_filename, 'w')
    f.write("id\tpbid\tis_fl\tstat\n")
    for line in open(group_filename):
        a,b=line.strip().split()
        for x in set(b.split(',')): f.write(f"{x}\t{a}\tY\tunique\n")
    f.close()
    print(f"Output written to: {output_filename}")

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("group_txt", help="Collapsed .group.txt to convert")
    parser.add_argument("output_txt", help="Output filename")

    args = parser.parse_args()
    if not os.path.exists(args.group_txt):
        print("Input .group.txt file {0} does not exist! Abort!".format(args.group_txt))
        sys.exit(-1)

    main(args.group_txt, args.output_txt)