import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tqdm import tqdm_notebook as tqdm
import re
from MulticoreTSNE import MulticoreTSNE as TSNE  # https://github.com/DmitryUlyanov/Multicore-TSNE


def plot_ts_like(df):
    """
simple plot with peaks in places where structure/nucl is
    :param df: dataframe with start and end columns
    :return: plt (for plt.show())
    """
    ranged = 100000
    cumsum_arr = np.zeros(ranged)
    for ind, row in df.iterrows():
        if int(row['end']) - 10000 > ranged:
            break
        for i in range(int(row['start']) - 10000, int(row['end']) - 10000):
            cumsum_arr[i] += 1
    x = np.arange(ranged)
    y = cumsum_arr

    plt.plot(x,y)
    # plt.show()
    return plt


def plot_tsne(ss10bp, position='struct'):
    """
form  physical features dataframe and use t-sne on it
https://lvdmaaten.github.io/tsne/
    :param position: before, after or struct for position relative to structure
    :param ss10bp: see preprocessing.py, sec. structs within 10bp from nucleosome
    :return: plt for plt.show()
    """
    diprodb = pd.read_csv('data/dprops.csv', index_col=0)

    features = []
    for struct in tqdm(ss10bp[position]):
        strl = re.findall('..', struct)
        temp = []
        for dyad in strl:
            temp.append(diprodb[diprodb['PropertyName'] == dyad].values.tolist()[0])
        features.append(pd.DataFrame(temp, columns=diprodb.columns).sum())
    features = pd.DataFrame(features)


    tsne = TSNE(n_jobs=8)
    transformed = tsne.fit_transform(features.loc[:, features.columns != 'PropertyName'])

    vis_x = transformed[:, 0]
    vis_y = transformed[:, 1]

    plt.scatter(vis_x, vis_y, cmap=plt.cm.get_cmap("jet", 10))
    plt.clim(-0.5, 9.5)
    return plt
    # plt.show()