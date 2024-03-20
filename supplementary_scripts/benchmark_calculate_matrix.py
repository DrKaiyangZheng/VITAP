import sys
import csv
from collections import defaultdict

'''input_file
# contig	Bona_fide_Phylum	Bona_fide_Class	Bona_fide_Order	Bona_fide_Family	Bona_fide_Genus	Predict_Phylum	Predict_Class	Predict_Order	Predict_Family	Predict_Genus
AB000906|_sliding:1-5000	Pisuviricota	Pisoniviricetes	Picornavirales	Iflaviridae	Iflavirus	Pisuviricota	Pisoniviricetes	Picornavirales	Iflaviridae	Iflavirus
AB002632|_sliding:1-5000	Hofneiviricota	Faserviricetes	Tubulavirales	Inoviridae	Saetivirus	NA	NA	NA	NA	NA
AB003289|_sliding:1-5000	Kitrinoviricota	Flasuviricetes	Amarillovirales	Flaviviridae	Pegivirus	Kitrinoviricota	Flasuviricetes	Amarillovirales	Flaviviridae	Pegivirus
AB003292|_sliding:1-5000	Kitrinoviricota	Flasuviricetes	Amarillovirales	Flaviviridae	Pegivirus	Kitrinoviricota	Flasuviricetes	Amarillovirales	Flaviviridae	Pegivirus
AB005247|_sliding:1-5000	Artverviricota	Revtraviricetes	Ortervirales	Metaviridae	Metavirus	Streptophyta	Magnoliopsida	Brassicales	Brassicaceae	Arabidopsis
AB005247|_sliding:10001-15000	Artverviricota	Revtraviricetes	Ortervirales	Metaviridae	Metavirus	Streptophyta	Magnoliopsida	Brassicales	Brassicaceae	Arabidopsis
AB005247|_sliding:15001-20000	Artverviricota	Revtraviricetes	Ortervirales	Metaviridae	Metavirus	Streptophyta	Magnoliopsida	Brassicales	Brassicaceae	Arabidopsis
AB005247|_sliding:20001-25000	Artverviricota	Revtraviricetes	Ortervirales	Metaviridae	Metavirus	Streptophyta	Magnoliopsida	Brassicales	Brassicaceae	Arabidopsis
AB005247|_sliding:25001-30000	Artverviricota	Revtraviricetes	Ortervirales	Metaviridae	Metavirus	Streptophyta	Magnoliopsida	Brassicales	Brassicaceae	-
'''

def compute_metrics(bona_fide, predict):
	TP, FN, FP, TN = 0, 0, 0, 0
	for bf, pr in zip(bona_fide, predict):
		if pr == "NA":  
			FN += 1
		elif bf == pr:  
			TP += 1
		else:  
			FP += 1

	return {"TP": TP, "FP": FP, "FN": FN, "TN": TN}

def main():
	input_file = sys.argv[1]
	output_file = sys.argv[2]

	phylum_aggregate = defaultdict(lambda: defaultdict(list))
	rank_aggregate = defaultdict(list)

	with open(input_file, 'r') as infile:
		reader = csv.DictReader(infile, delimiter="\t")
		for row in reader:
			bona_fide_phylum = row['Bona_fide_Phylum']
			if bona_fide_phylum != '-':
				for rank in ['Phylum', 'Class', 'Order', 'Family', 'Genus']:
					bona_fide = row[f'Bona_fide_{rank}']
					predict = row[f'Predict_{rank}']

					metrics = compute_metrics([bona_fide], [predict])

					phylum_aggregate[bona_fide_phylum][rank].append(metrics)
					rank_aggregate[rank].append(metrics)

	# Aggregate metrics
	def aggregate_metrics(metrics_list):
		TP = sum(m["TP"] for m in metrics_list)
		FP = sum(m["FP"] for m in metrics_list)
		FN = sum(m["FN"] for m in metrics_list)
		TN = sum(m["TN"] for m in metrics_list)
		
		Accuracy = (TP + TN) / (TP + TN + FP + FN) if TP + TN + FP + FN != 0 else 0
		Precision = TP / (TP + FP) if TP + FP != 0 else 0
		Recall = TP / (TP + FN) if TP + FN != 0 else 0
		F1_score = 2 * (Precision * Recall) / (Precision + Recall) if Precision + Recall != 0 else 0
		
		return {"Accuracy": Accuracy, "Precision": Precision, "Recall": Recall, "F1_score": F1_score}

	# Write to files
	phylum_output_file = output_file + "_phylum.tsv"
	rank_output_file = output_file + "_rank.tsv"

	with open(phylum_output_file, 'w') as pout:
		writer = csv.DictWriter(pout, fieldnames=['Phylum', 'Rank', 'Accuracy', 'Precision', 'Recall', 'F1_score'], delimiter='\t')
		writer.writeheader()
		for phylum, ranks in phylum_aggregate.items():
			for rank, metrics_list in ranks.items():
				aggregated = aggregate_metrics(metrics_list)
				writer.writerow({**{"Phylum": phylum, "Rank": rank}, **aggregated})

	with open(rank_output_file, 'w') as rout:
		writer = csv.DictWriter(rout, fieldnames=['Rank', 'Accuracy', 'Precision', 'Recall', 'F1_score'], delimiter='\t')
		writer.writeheader()
		for rank, metrics_list in rank_aggregate.items():
			aggregated = aggregate_metrics(metrics_list)
			writer.writerow({**{"Rank": rank}, **aggregated})

if __name__ == '__main__':
	main()
