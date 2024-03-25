import argparse

def createplotstats(coverage):
    '''
    Generate necessary stats to plot CANOES analysis bar and boxplot
    Return a dictionary with stats
    '''
    dictionary = {'header': [], 'depth': []}
    with open(coverage, 'r') as file:
        for line in file:
            if line.startswith('chrom'):
                header_line = line.strip().split('\t')
                dictionary['header'].append('Chr:Start-End')
                dictionary['header'].extend(header_line[3:])
            elif line.startswith('chrX') or line.startswith('chrY'):
                continue
            else:
                line_elements = line.strip().split('\t')
                metrics = []
                for depth in map(int, line_elements[3:]):
                    try:
                        metrics.append(str(round(float(depth / (int(line_elements[2]) - int(line_elements[1]))), 4)).replace(',', '.'))
                    except ZeroDivisionError:
                        metrics.append('0')
                dictionary['depth'].append(line_elements[0] + ':' + line_elements[1] + '-' + line_elements[2] + '\t' + '\t'.join(metrics))
    return dictionary

def main():
    parser = argparse.ArgumentParser(description='Generate statistics for CANOES analysis')
    parser.add_argument('input', help='Input coverage file')
    parser.add_argument('output', help='Output file to write the stats')
    args = parser.parse_args()

    values = createplotstats(args.input)
    with open(args.output, 'w+') as o:
        o.write('\t'.join(values['header'])+'\n')
        for val in values['depth']:
            o.write(val+'\n')
    print("#[INFO] Ending parse "+args.input+" write output in "+args.output)

if __name__ == "__main__":
    main()