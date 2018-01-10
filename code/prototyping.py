# get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os
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

    df = pd.concat([df1, df2])
    df.dropna(inplace=True)
    df.to_csv('chr1.S15-30.filtered.csv')
    return df


def s15_30(uplifted, df):
    """
make Nucleosome-centered positions of sec. structures plot, counts positions of sec.structures' nucleotides
    """

    cumsum_arr = np.zeros(1000)

    for peak_pos in tqdm(pd.to_numeric(uplifted.peak_pos)):
        temp_df = df[(df['start'] > peak_pos - 500) & (df['end'] < peak_pos + 500)]
        for ind, row in temp_df.iterrows():
            for i in range(int(row['start']) - peak_pos + 500, int(row['end']) - peak_pos + 500):
                cumsum_arr[i] += 1

    x = np.arange(-500, 500)
    plt.plot(x, cumsum_arr)
    plt.ylabel('Number of sec. structures')
    plt.xlabel('Nucleosome-centered positions of sec. structures')
    plt.savefig('Chr1_S15-30_centered_plot.png')


# ### Because S15-30 gave irrelevant result, we'll test S16-50


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
    pd.DataFrame(arr_4_rows).to_csv('sec_struct_0-10bp_to_nucleosome.csv', index=False)
