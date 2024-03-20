import pandas as pd
import sys

'''
#input_file
UViGs	NR_Hit_Family	NR_Hit_Family_Count	VITAP	IMGVR
IMGVR_UViG_2029527005_000017	Peduoviridae	0.857142857	Peduoviridae	-
IMGVR_UViG_2037747001_000002	Anaerodiviridae	0.416666667	Anaerodiviridae	-
IMGVR_UViG_2037747001_000002	Haloferuviridae	0.027777778	Anaerodiviridae	-
IMGVR_UViG_2043231000_000017	Corticoviridae	0.25	Peduoviridae	-
IMGVR_UViG_2043231000_000017	Peduoviridae	0.5	Peduoviridae	-
IMGVR_UViG_2049941000_000109	Ackermannviridae	0.333333333	Kyanoviridae	-
IMGVR_UViG_2051223000_000015	Intestiviridae	0.857142857	Intestiviridae	-
IMGVR_UViG_2081372001_000001	Peduoviridae	0.833333333	Peduoviridae	-
...
#output_file
UViGs	NR_Hit_Family	NR_Hit_Family_Count	VITAP	IMGVR	NR_Top_Family	VITAP_agreement	IMGVR_agreement
IMGVR_UViG_2029527005_000017	Peduoviridae	0.857142857	Peduoviridae	-	Peduoviridae	√	
IMGVR_UViG_2037747001_000002	Anaerodiviridae	0.416666667	Anaerodiviridae	-	Anaerodiviridae	√	
IMGVR_UViG_2043231000_000017	Corticoviridae	0.25	Peduoviridae	-	Peduoviridae	√	
IMGVR_UViG_2049941000_000109	Ackermannviridae	0.333333333	Kyanoviridae	-	Ackermannviridae		
IMGVR_UViG_2051223000_000015	Intestiviridae	0.857142857	Intestiviridae	-	Intestiviridae	√	
IMGVR_UViG_2081372001_000001	Peduoviridae	0.833333333	Peduoviridae	-	Peduoviridae	√	
IMGVR_UViG_2084038008_000017	Autographiviridae	0.111111111	Rhodogtaviriformidae	-	Autographiviridae;Casjensviridae;Mesyanzhinovviridae		
IMGVR_UViG_2084038019_000043	Druskaviridae	0.5	Druskaviridae	-	Druskaviridae	√	
IMGVR_UViG_2088090007_000003	Circoviridae	0.666666667	Circoviridae	-	Circoviridae	√	
'''

def process_tsv(input_file, output_file):
    # loading TSV file
    df = pd.read_csv(input_file, sep='\t')

    # For each UViG, find all the NR_Hit_Families corresponding to the maximum value of NR_Hit_Family_Count
    def get_top_families(group):
        max_count = group['NR_Hit_Family_Count'].max()
        return ';'.join(group[group['NR_Hit_Family_Count'] == max_count]['NR_Hit_Family'].unique())

    top_families = df.groupby('UViGs').apply(get_top_families)
    df = df.merge(top_families.rename('NR_Top_Family'), on='UViGs')

    # Check the VITAP and IMGVR columns
    df['VITAP_agreement'] = df.apply(lambda x: '√' if x['VITAP'] in x['NR_Top_Family'].split(';') else '', axis=1)
    df['IMGVR_agreement'] = df.apply(lambda x: '√' if x['IMGVR'] in x['NR_Top_Family'].split(';') else '', axis=1)

    # Saving processed result
    df.drop_duplicates(subset='UViGs').to_csv(output_file, sep='\t', index=False)

if __name__ == "__main__":
    process_tsv(sys.argv[1], sys.argv[2])
