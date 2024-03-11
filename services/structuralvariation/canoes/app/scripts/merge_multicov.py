import argparse
import pandas as pd

def generate_multicov_output(input_files, output_file):
    # Initialize an empty DataFrame
    df = pd.DataFrame()

    # Define the header
    header = ["chrom", "start", "end"]

    # Loop through input files
    for multicovfile in input_files:
        sample_name = multicovfile.split("/")[-1].split(".")[0]
        header.append(sample_name)

        # Read data from file and append to DataFrame
        with open(multicovfile, "r") as file:
            for line in file:
                fields = line.strip().split("\t")
                index = "\t".join(fields[:3])
                value = fields[-1]

                # Check if index already exists in DataFrame
                if index not in df.index:
                    df.loc[index] = [None] * len(header)
                df.at[index, sample_name] = value

    # Write DataFrame to output file
    df.columns = header
    df.to_csv(output_file, sep='\t', index=False)

def main():
    parser = argparse.ArgumentParser(description='Generate multicov output')
    parser.add_argument('input_files', nargs='+', help='Input multicov files')
    parser.add_argument('-o', '--output', help='Output file')
    args = parser.parse_args()

    if not args.input_files:
        print("Error: At least one input file is required.")
        return

    if not args.output:
        print("Error: Output file path is required.")
        return

    generate_multicov_output(args.input_files, args.output)

if __name__ == "__main__":
    main()