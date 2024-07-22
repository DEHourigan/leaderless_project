import sys

def process_files(file_paths, evalue_threshold):
    results = {}
    for file_path in file_paths:
        with open(file_path, 'r') as file:
            for line in file:
                if line.startswith('#'):
                    continue
                parts = line.split()
                try:
                    evalue = float(parts[4])
                except IndexError:
                    continue
                if evalue < evalue_threshold:
                    seq_id = parts[0]
                    if seq_id not in results or results[seq_id][4] > evalue:
                        results[seq_id] = parts
    sorted_results = sorted(results.values(), key=lambda x: (x[0], float(x[4])))
    for result in sorted_results:
        print("\t".join(result))

if __name__ == "__main__":
    file_paths = sys.argv[1:-1]
    evalue_threshold = float(sys.argv[-1])
    process_files(file_paths, evalue_threshold)

