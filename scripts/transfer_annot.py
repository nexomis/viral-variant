from Bio import SeqIO
from Bio import Align
import argparse
import sys
import csv

def get_ref2alt(ref_file, alt_file, base1, verbose):
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 5
    aligner.mismatch_score = -4
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    aligner.target_end_gap_score = -5
    aligner.query_end_gap_score = -5
    parsed_records = dict()
    parsed_records["ref"] = list(SeqIO.parse(ref_file, "fasta"))
    parsed_records["alt"] = list(SeqIO.parse(alt_file, "fasta"))
    records = dict()
    for key in parsed_records.keys():
        if len(parsed_records[key]) != 1:
            sys.exit("1 and only 1 contig per fasta file required")
        records[key] = parsed_records[key][0]
    if verbose:
        print("aligning %s with %s" % (records["ref"].id, records["alt"].id))
    alignments = aligner.align(records["ref"].seq, records["alt"].seq)
    best = alignments[0]
    if verbose:
        # Print the header for the summary table
        print("Alignment Length | Aligned Bases from Target | Aligned Bases from Query | Alignment Score")
        print("-" * 80)
        # Iterate through the alignments and print the details
        for alignment in alignments:
            alignment_array = alignment.__array__()
            target = "".join(map(lambda x: x.decode("utf-8"), alignment_array[0]))
            query = "".join(map(lambda x: x.decode("utf-8"), alignment_array[1]))
            alignment_length = len(target)
            target_aligned_bases = sum(1 for x in target if x != '-')
            query_aligned_bases = sum(1 for x in query if x != '-')
            alignment_score = alignment.score

            print(f"{alignment_length:^18} | {target_aligned_bases:^25} | {query_aligned_bases:^27} | {alignment_score:^18}")
    alignment_array = best.__array__()

    ref = "".join(map(lambda x: x.decode("utf-8"), alignment_array[0]))
    alt = "".join(map(lambda x: x.decode("utf-8"), alignment_array[1]))
    middle = ""
    for char1, char2 in zip(alignment_array[0], alignment_array[1]):
        if char1 == char2 and char1 != b'-':
            middle += '|'
        elif char1 == b'-' or char2 == b'-':
            middle += ' '
        else:
            middle += '.'
    # start at 0 so that if no gap then the start seq is 1


    if base1:
        pos_alt = 0
        ref2alt = [pos_alt]
    else:
        pos_alt = -1
        ref2alt = []
    
    pos_ref_al = -1 # for verbose only
    pos_alt_al = -1 # for verbose only

    al2ref = []
    al2alt = []

    for i in range(0, len(ref)):
        if alt[i] != "-":
            pos_alt += 1
            pos_alt_al += 1
        if ref[i] != "-":
            pos_ref_al += 1
            ref2alt.append(pos_alt)
        al2ref.append(pos_ref_al)
        al2alt.append(pos_alt_al)
    # add another iteration for non inclusive end
    ref2alt.append(pos_alt + 1)

    if verbose:
        win_size = 50 # print window size
        i = 0
        s = win_size * (i)
        e = win_size * (i + 1)

        def format_pos(pos):
            return str(pos) + " " * (10 - len(str(pos)))

        while ref[s:e] != "":
            s = win_size * (i)
            e = win_size * (i + 1)
            if e > len(ref):
                e = len(ref)
            mid_s = format_pos(s)
            mid_e = " " + str(e - 1)
            ref_s = format_pos(al2ref[s])
            ref_e = " " + str(al2ref[e - 1])
            alt_s = format_pos(al2alt[s])
            alt_e = " " + str(al2alt[e - 1])
            print("")
            print("print alignment with 0-based coord")
            print("")
            print(ref_s + ref[s:e] + ref_e)
            print(mid_s + middle[s:e] + mid_e)
            print(alt_s + alt[s:e] + alt_e)
            i = i + 1
    return ref2alt

def transform_coordinates(source_index, coordinate_mapping):
    len_map = len(coordinate_mapping)
    if source_index >= len_map:
        raise ValueError(f"Index ({source_index}) >= List size {len_map}")
    return max(coordinate_mapping[source_index], 0)

def process_file(input_file, output_file, columns_to_transform,    start_column,
    size_column, column_delimiter, intra_column_delimiter, coordinate_mapping):
    with open(input_file, "r") as infile, \
        open(output_file, "w", newline="") as outfile:
        reader = csv.reader(infile, delimiter=column_delimiter)
        writer = csv.writer(outfile, delimiter=column_delimiter)

        for row in reader:
            for col in columns_to_transform:
                if col >= len(row):
                    raise ValueError("Given column not present")
                values = row[col].split(intra_column_delimiter)
                transformed_values = [transform_coordinates(int(value), coordinate_mapping) for value in values]
                row[col] = intra_column_delimiter.join(map(str, transformed_values))

            if start_column is not None and size_column is not None and start_column < len(row) and size_column < len(row):
                start_values = row[start_column].split(intra_column_delimiter)
                size_values = row[size_column].split(intra_column_delimiter)
                transformed_start_values = []
                transformed_size_values = []

                for start, size in zip(start_values, size_values):
                    start_value = int(start)
                    size_value = int(size)
                    end_value = start_value + size_value
                    transformed_start = transform_coordinates(start_value, coordinate_mapping)
                    transformed_end = transform_coordinates(end_value, coordinate_mapping)
                    transformed_size = transformed_end - transformed_start
                    transformed_start_values.append(str(transformed_start))
                    transformed_size_values.append(str(transformed_size))

                row[start_column] = intra_column_delimiter.join(transformed_start_values)
                row[size_column] = intra_column_delimiter.join(transformed_size_values)

            writer.writerow(row)

def main():
    parser = argparse.ArgumentParser(description = "Transfer an annotation " + \
        "from one contig to another (in whatever format)")
    parser.add_argument("-ref", metavar = "fasta", help = "reference fasta")
    parser.add_argument("-alt", metavar = "fasta", help = "alternative fasta")
    parser.add_argument("-annot", metavar = "gff",
        help = "annotation file in GFF format")
    parser.add_argument("-out", metavar = "gff",
        help = "annotation for alternative fasta")
    parser.add_argument("--verbose", action = "store_true",
        help = "Enable verbose logging")
    parser.add_argument("--base1", action = "store_true",
        help = "First base is 1 not 0")
    parser.add_argument("--columns", default = "1,2,6,7",
        help = "Columns to be replaced, e.g., \"2,5,9\".")
    parser.add_argument("--start_column", default = 11, type = int, 
        help = "Column representing the start coordinate.")
    parser.add_argument("--size_column", default = 10, type = int,
        help = "Column representing the size, related to the start column.")
    parser.add_argument("--column_delimiter", default = "\t",
        help = "Delimiter used between columns.")
    parser.add_argument("--intra_column_delimiter", default = ",",
        help ="Delimiter used within columns.")

    args = parser.parse_args()

    ref2alt = get_ref2alt(args.ref, args.alt, args.base1, args.verbose)

    columns_to_transform = list(map(int, args.columns.split(",")))

    process_file(args.annot, args.out, columns_to_transform, args.start_column, args.size_column, args.column_delimiter, args.intra_column_delimiter, ref2alt)

if __name__ == "__main__":
    main()


