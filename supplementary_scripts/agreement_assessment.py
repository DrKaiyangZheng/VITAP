import sys

'''
file A and B
IMGVR_UViG_2529293226_000001	Faserviricetes
IMGVR_UViG_2504756089_000001	Faserviricetes
IMGVR_UViG_2505313053_000001	Faserviricetes
IMGVR_UViG_2506210030_000001	Faserviricetes
IMGVR_UViG_2767802255_000008	Faserviricetes
IMGVR_UViG_2698536412_000001	Faserviricetes
IMGVR_UViG_2510461032_000002	Faserviricetes
IMGVR_UViG_2516653000_000014	Faserviricetes
IMGVR_UViG_2516653000_000002	Faserviricetes
IMGVR_UViG_2517572009_000001	-
...

The newly identified items are based on file B. The '-' value in file B were exclued in agreement assessment

Use the first column as the index (genome), to check whether the values in the second column are consistent across two files. If there are duplicate items in the first column, check whether all corresponding values in the second column for these duplicates are consistent across the two files. According to this rule, count the number of genomes for which the second column values remain consistent in both files:
1) The consistency check of the second column's values only applies to genomes in file B where the second column value is not "-". If the second column value for a genome in file B is "-", then it is not subjected to consistency verification and is not included in the final count;
2)For genomes in file B with a second column value of "-", if the corresponding genome in file A has a second column value that is not "-", then include it in the "total number of newly identified items" count, and return it as a third type of result to the terminal;
3) If a genome in both file A and file B has the second column value as '-', then also include it in the consistency count;
4) For genomes in file A with a second column value of "-", if the corresponding genome in file B has a second column value that is not "-", then include it in the "total number of lost identified items" count, and return it as a fourth type of result to the terminal.
'''

def read_tsv_to_dict(filename):
    data_dict = {}
    with open(filename, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            genome = parts[0]
            value = parts[1] if len(parts) > 1 else "-"
            if genome not in data_dict:
                data_dict[genome] = []
            data_dict[genome].append(value)
    return data_dict

def compare_dicts(dictA, dictB):
    agreement_count = 0
    newly_identified_count = 0
    lost_identified_count = 0
    disagreement_count = 0
    discrepancies = []
    newly_identified_items = []

    total_count = len(set(dictA.keys()).union(set(dictB.keys())))

    for genome, valuesB in dictB.items():
        if genome not in dictA:
            continue

        valuesA = dictA[genome]

        if all(v == "-" for v in valuesA) and any(v != "-" for v in valuesB):
            lost_identified_count += 1
        elif all(v == "-" for v in valuesB) and all(v == "-" for v in valuesA):
            agreement_count += 1
        elif all(v == "-" for v in valuesB) and any(v != "-" for v in valuesA):
            newly_identified_count += 1
            newly_identified_items.append((genome, ', '.join(valuesA), ', '.join(valuesB)))
        elif any(v != "-" for v in valuesA) and any(v != "-" for v in valuesB):
            if set(valuesA) == set(valuesB):
                agreement_count += 1
            else:
                disagreement_count += 1
                discrepancies.append((genome, ', '.join(valuesA), ', '.join(valuesB)))

    return agreement_count, total_count, newly_identified_count, lost_identified_count, disagreement_count, discrepancies, newly_identified_items

if __name__ == '__main__':
    if len(sys.argv) != 5:
        print("Usage: python script_name.py fileA.tsv fileB.tsv discrepancies.txt newly_identified.txt")
        sys.exit(1)

    fileA = sys.argv[1]
    fileB = sys.argv[2]
    discrepancies_file = sys.argv[3]
    newly_identified_file = sys.argv[4]

    dictA = read_tsv_to_dict(fileA)
    dictB = read_tsv_to_dict(fileB)

    agreement, total, newly_identified, lost_identified, disagreement, discrepancies_list, newly_identified_list = compare_dicts(dictA, dictB)

    with open(discrepancies_file, 'w') as file:
        for genome, valueA, valueB in discrepancies_list:
            file.write(f"{genome}\t{valueA}\t{valueB}\n")

    with open(newly_identified_file, 'w') as file:
        for genome, valueA, valueB in newly_identified_list:
            file.write(f"{genome}\t{valueA}\t{valueB}\n")

    print(f"total number of items: {total}")
    print(f"number of items with agreement: {agreement}") 
    print(f"number of items with disagreement: {disagreement}")     
    print(f"total number of newly identified items: {newly_identified}")
    print(f"total number of lost identified items: {lost_identified}")










