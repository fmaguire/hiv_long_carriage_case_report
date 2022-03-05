#!/usr/bin/env python
import re
import sys
import pandas as pd

if __name__ == "__main__":

    dnds_data = {'Seq1': [], 'Seq2': [], 'Ka/Ks': []}
    with open(sys.argv[1], "r") as fh:
        for line in fh:
            line = line.strip()
            if len(line) == 0:
                continue

            if re.match("^\d+ \(.*", line):
                dnds_data['Seq1'].append(line.split('(')[1].split(')')[0])
                dnds_data['Seq2'].append(line.split('(')[2].split(')')[0])

                next(fh)
                next(fh)
                next(fh)
                dnds_line = next(fh)
                dnds_data['Ka/Ks'].append(float(dnds_line.strip().split()[7]))
    dnds_data = pd.DataFrame(dnds_data)
    dnds_data.to_csv('codeml_kaks_values.tsv', sep='\t', index=False)
