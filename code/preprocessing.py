# get_ipython().run_line_magic('matplotlib', 'inline')
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tqdm import tqdm_notebook as tqdm


def make_uplift_csv():
    # Preparing data from web
    with open('peaks.min_peak_score_0.6.thresh_0.5.txt', 'r') as f:
        peaks = f.read()
    # [4:] is because of the annotation lines
    df_peaks = pd.DataFrame([sub.split("\t") for sub in peaks.split('\n')][4:],
                            columns=['chrom', 'start', 'end', 'peak_pos', 'score'])
    df_peaks.drop(df_peaks.tail(1).index, inplace=True)
    df_peaks.head()

    # https://pypi.python.org/pypi/pyliftover
    from pyliftover import LiftOver
    lo = LiftOver('hg18', 'hg19')
    lo.convert_coordinate('chr1', 1000000)  # testing 

    row_to_delete = []

    def liftover(row):
        value = lo.convert_coordinate(row['chrom'], int(row['start']))
        if not value:
            row_to_delete.append(row)
        else:
            row['start'] = int(value[0][1])

        value = lo.convert_coordinate(row['chrom'], int(row['end']))
        if not value:
            row_to_delete.append(row)
        else:
            row['end'] = int(value[0][1])

        value = lo.convert_coordinate(row['chrom'], int(row['peak_pos']))
        if not value:
            row_to_delete.append(row)
        else:
            row['peak_pos'] = int(value[0][1])

        return row

    uplifted = df_peaks.apply(liftover, axis=1)

    # Deleting rows that can't be uplifted
    ind_row_to_delete = [x.name for x in row_to_delete]
    ind_row_to_delete = list(set(ind_row_to_delete))
    len(ind_row_to_delete)
    uplifted.drop(ind_row_to_delete, inplace=True)
    uplifted.to_csv('peaks.min_peak_score_0.6.thresh_0.5.csv')
    print(uplifted.tail())
    return uplifted


def make_filtered_s15_30():
    """
make clean dataframe with 5 columns almost without duplicates
Unnamed: 0	    start	end	      len_stem	       len_loop
    :return: this df
    """
    base_path = "/home/shared/STEMLOOPS/hg19/S15-30_L0-10_M5"
    filename = "chr1.fna.S15-30_L0-10_M5.pal"
    path_to_file = os.path.join(base_path, filename)
    with open(path_to_file, 'r') as f:
        #     next(f)  # if there is a description line
        splitted = f.read().split('\n')  # raw file separated by \n's

    from tqdm import tqdm_notebook as tqdm
    cols = ['start', 'end', 'len_stem', 'len_loop']  # , 'seq_left', 'seq_right', 'full_seq', '1', '2', '3']
    # will do in 2 parts because file is too big

    df1 = pd.DataFrame([sub.split("\t")[:4] for sub in tqdm(splitted[:10000000])], columns=cols)
    df1 = df1.apply(pd.to_numeric)
    df1 = df1.drop_duplicates(subset=['start']).drop_duplicates(['end'])

    df2 = pd.DataFrame([sub.split("\t")[:4] for sub in tqdm(splitted[10000000:])], columns=cols)
    df2 = df2.apply(pd.to_numeric)
    df2 = df2.drop_duplicates(subset=['start']).drop_duplicates(['end'])

    df = pd.concat([df1, df2]).reset_index(drop=True)
    df.dropna(inplace=True)
    df.to_csv('chr1.S15-30.filtered.csv')
    return df


def s16_50(uplifted):
    """
see s15_30
second part - Structure-centered positions of nucleosomes, counts positions of nucleosome's nucleotides
    """
    base_path = "/home/shared/STEMLOOPS/hg19/S16-50_L0-10_M3"
    filename = "chr1.fa.S16-50_L0-10_M3.pal.cleaned"
    path_to_file = os.path.join(base_path, filename)
    with open(path_to_file, 'r') as f:
        #     next(f)  # if there is description line
        splitted = f.read().split('\n')  # raw file separated by \n's

    cols = ['start', 'end', 'len_stem', 'len_loop']  # , 'seq_left', 'seq_right', 'full_seq', '1', '2', '3']
    df = pd.DataFrame([sub.split("\t")[:4] for sub in tqdm(splitted)], columns=cols)
    df = df.apply(pd.to_numeric)
    df.drop(df.tail(1).index, inplace=True)
    df['center'] = df[['start', 'end']].mean(axis=1).astype('int')
    df.tail()

    cumsum_arr = np.zeros(1000)

    for peak_pos in tqdm(pd.to_numeric(uplifted.peak_pos)):
        temp_df = df[(df['start'] > peak_pos - 500) & (df['end'] < peak_pos + 500)]
        for ind, row in temp_df.iterrows():
            for i in range(int(row['start']) - peak_pos + 500, int(row['end']) - peak_pos + 500):
                cumsum_arr[i] += 1
    # cumsum_arr

    x = np.arange(-500, 500)
    plt.plot(x, cumsum_arr)
    plt.ylabel('Number of sec. structures')
    plt.xlabel('Nucleosome-Centered position of sec. structures')
    plt.savefig('pictures/Chr1_S16-50_centered_plot.png')

    # Let's reverse and center sec.structures and count nucleosomes around
    uplifted = uplifted[uplifted['chrom'] == 'chr1']
    uplifted.drop('chrom', axis=1, inplace=True)
    uplifted = uplifted.apply(pd.to_numeric)
    uplifted.tail()

    cumsum_arr = np.zeros(1000)

    for peak_pos in tqdm(df.center):
        temp_df = uplifted[(uplifted['start'] > peak_pos - 500) & (uplifted['end'] < peak_pos + 500)]
        for ind, row in temp_df.iterrows():
            for i in range(int(row['start']) - peak_pos + 500, int(row['end']) - peak_pos + 500):
                cumsum_arr[i] += 1
    # cumsum_arr

    x = np.arange(-500, 500)
    plt.plot(x, cumsum_arr)
    plt.ylabel('Number of nucleosomes')
    plt.xlabel('Structure(S16-50)-centered_nucleosome_positions')
    plt.savefig('pictures/structure(S16-50)-centered_plot_Chr1_nucleosome_positions.png')

# ### It's better. Let's check quadruplexes:


def quadr(uplifted):
    """
see s15_30
    :param uplifted: uplifted array from function above
    """
    base_path = "/home/konovalovdmitry/results"
    filename = "chr1.out"
    path_to_file = os.path.join(base_path, filename)
    df = pd.read_csv(path_to_file, sep=' ', header=None)
    cols = ['chr', 'start', 'end', 'seq']
    df.columns = cols
    df.drop(['seq', 'chr'], axis=1, inplace=True)
    df = df.apply(pd.to_numeric)

    # # or you can open it like this:
    # all_quad_coords = pd.read_csv('/home/konovalovdmitry/results/all_quad_coords.tsv', sep='\t', header=None)
    # all_quad_coords = all_quad_coords[all_quad_coords[0] == 'chr1']
    # cols = ['chr', 'start', 'end']
    # all_quad_coords.columns = cols
    # df = all_quad_coords.copy()

    cumsum_arr = np.zeros(1000)

    for peak_pos in tqdm(pd.to_numeric(uplifted.peak_pos)):
        temp_df = df[(df['start'] > peak_pos - 500) & (df['end'] < peak_pos + 500)]
        for ind, row in temp_df.iterrows():
            for i in range(int(row['start']) - peak_pos + 500, int(row['end']) - peak_pos + 500):
                cumsum_arr[i] += 1
    # cumsum_arr

    x = np.arange(-500, 500)
    plt.plot(x, cumsum_arr)
    plt.ylabel('Number of sec. structures')
    plt.xlabel('Centered position')
    plt.savefig('Chr1_quadruplexes_centered_plot.png')


def make_10bp(df, uplifted):
    """
Выделить все структуры, у которых основание ножки находятся на расстоянии от 0 до 10 нуклеотидов от границы нуклеосомы.
    :param df: sec. struct. df
    :param uplifted: nucleosome df
    """
    arr_4_rows = []
    for start_upl in tqdm(pd.to_numeric(uplifted.start)):
        temp_df = df[(df['end'] > start_upl - 10) & (df['end'] < start_upl)]
        if not temp_df.empty:
            for ind, row in temp_df.iterrows():
                arr_4_rows.append(row)
    ss10bp = pd.DataFrame(arr_4_rows)
    ss10bp.to_csv('sec_struct_0-10bp_to_nucleosome.csv', index=False)
    return ss10bp.reset_index(drop=True)


def make_n_bp(df, uplifted, n=70):
    """
Return all 'df' elements where end in range n before start of 'uplifted'
    :param n: max number of bp from df.end to uplifted.start
    :param df: sec. struct. df with columns start   end len_stem    len_loop (10941426 rows)
    :param uplifted: nucleosome df with columns chrom   start   end peak_pos    score (1037956 rows)
    """
    df1 = df.sort_values('end').reset_index(drop=True)
    # getting indexes where should df.end be, if inserted into uplifted
    indexes = np.searchsorted(uplifted.start.values, df1.end.values)
    uplifted1 = uplifted.append(uplifted.tail(1))  # for df.end > uplifted.start.max()
    uplifted_ind = uplifted1.iloc[indexes].reset_index(drop=True)  # getting all uplifted.start closest to df.end
    uplifted_ind['df_end'] = df1.end.values
    uplifted_ind['diff'] = uplifted_ind.start - uplifted_ind.df_end
    mask = (uplifted_ind['diff'] >= 0) & (uplifted_ind['diff'] < n)
    return df1[mask].reset_index(drop=True)


def make_go_terms(ss10bp, suffix='quad'):
    """
makes go_terms.csv for pasting to http://revigo.irb.hr/ and getting genes and their types
also makes relevant_goa_names.csv for names of this genes
    """
    ptt = pd.read_csv('../data/ptt_hg19.txt', delimiter='\t')
    ptt1 = ptt[ptt['chrom'] == 'chr1']

    def make_10bp_pr(df, uplifted):
        arr_4_rows = []
        for start_upl in tqdm(pd.to_numeric(uplifted.end)):
            temp_df = df[(df['txStart'] - 1000 < start_upl) & (df['txEnd'] > start_upl) & (df.strand == '+')]
            if not temp_df.empty:
                for ind, row in temp_df.iterrows():
                    arr_4_rows.append(row)
        return pd.DataFrame(arr_4_rows)

    relevant_ptt = make_10bp_pr(ptt1, ss10bp)

    goa = pd.read_csv('../data/goa_human.gaf', delimiter='\t', header=None)

    relevant_goa = goa[goa[1].isin(relevant_ptt.proteinID)].drop_duplicates(1)
    relevant_goa[9].to_csv('../data/relevant_goa_names_{}.csv'.format(suffix))
    relevant_goa[4].to_csv('../data/go_terms_{}.csv'.format(suffix), index=False)


def add_bp_according2_start_end(neg_example):
    """
adds nucleotides string according to start-end columns
    :param neg_example: DataFrame with start-end columns
    """
    from Bio import SeqIO
    first_record = str(next(SeqIO.parse("/home/shared/hg19/chr1.fna", "fasta")).seq)

    temp = pd.DataFrame(np.ndarray((neg_example.shape[0], 4)),
                        columns=["struct", "before", "after", 'half_struct'], dtype='str')
    for idx, row in neg_example.iterrows():
        temp.loc[idx][0] = first_record[row['start']:row['end']]
        temp.loc[idx][1] = first_record[row['start'] - 20:row['start']]
        temp.loc[idx][2] = first_record[row['end']:row['end']    + 20]
        temp.loc[idx][3] = first_record[row['start']:row['start'] + (row['end'] - row['start'])//2]

    return pd.concat([neg_example, temp], axis=1)


def make_neg_example_and_concat(df, quad10with_bps):
    """
    :param df: DataFrame with start-end coordinates
    :param quad10with_bps: second DataFrame to be concated (ie ss10bp)
    :return: concated df for feature calculations
    """
    neg_example = df.sample(5000).reset_index(drop=True)
    neg_example = add_bp_according2_start_end(neg_example)
    neg_example['value'] = 0
    quad10with_bps['value'] = 1
    concated = pd.concat([neg_example, quad10with_bps])
    return concated


def calc_feats(concated, position='struct'):
    """
calculate features from diprodb
    :param position: struct, before, after or half_struct
    :param concated: df with bp as after add_bp_according2_start_end
    """
    diprodb = pd.read_csv('../data/dprops.csv', index_col=0)
    import re
    features = []
    for struct in tqdm(concated.loc[:, concated.columns != 'value'][position]):
        strl = re.findall('..', struct)
        temp = []
        for dyad in strl:
            temp.append(diprodb[diprodb['PropertyName'] == dyad].values.tolist()[0])
        features.append(pd.DataFrame(temp, columns=diprodb.columns).sum())
    return pd.DataFrame(features)
