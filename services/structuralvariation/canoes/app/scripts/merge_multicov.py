import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description="Process multicov files.")
    parser.add_argument("input", nargs="+", help="Input multicov files")
    parser.add_argument("-o", "--output", required=True, help="Output file")
    return parser.parse_args()

def main():
    args = parse_arguments()
    
    data = {}
    header = ["chrom", "start", "end"]
    
    for multicovfile in args.input:
        sample_name = multicovfile.split("/")[-1].split(".")[2]
        header.append(sample_name)
        
        with open(multicovfile, "r") as i:
            for line in i:
                field = line.strip().split("\t")
                index = "\t".join(field[0:3])
                value = field[-1]
                
                if index not in data:
                    data[index] = {}
                data[index][sample_name] = value

    with open(args.output, "w") as o:
        o.write("\t".join(header) + "\n")
        for index in data:
            sample_value = [index]
            for sample_name in header[3:]:
                sample_value.append(data[index].get(sample_name, "NA"))
            o.write("\t".join(sample_value) + "\n")

if __name__ == "__main__":
    main()