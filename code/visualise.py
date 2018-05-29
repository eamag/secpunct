import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import six
from MulticoreTSNE import MulticoreTSNE as TSNE  # https://github.com/DmitryUlyanov/Multicore-TSNE
from tqdm import tqdm_notebook as tqdm


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


def plot_tsne(features, position='struct'):
    """
form  physical features dataframe and use t-sne on it
https://lvdmaaten.github.io/tsne/
    :param features: DataFrame of features
    :param position: before, after or struct for position relative to structure
    :param ss10bp: (deprecated) see preprocessing.py, sec. structs within 10bp from nucleosome
    :return: plt for plt.show()
    """
    diprodb = pd.read_csv('../data/dprops.csv', index_col=0)

    # features = []
    # for struct in tqdm(ss10bp[position]):
    #     strl = re.findall('..', struct)
    #     temp = []
    #     for dyad in strl:
    #         temp.append(diprodb[diprodb['PropertyName'] == dyad].values.tolist()[0])
    #     features.append(pd.DataFrame(temp, columns=diprodb.columns).sum())
    # features = pd.DataFrame(features)


    tsne = TSNE(n_jobs=8)
    transformed = tsne.fit_transform(features.loc[:, features.columns != 'PropertyName'])

    vis_x = transformed[:, 0]
    vis_y = transformed[:, 1]

    plt.scatter(vis_x, vis_y, cmap=plt.cm.get_cmap("jet", 10))
    plt.clim(-0.5, 9.5)
    return plt
    # plt.show()
    # or we can do this and label relevant vs random
    # diprodb = pd.read_csv('../data/dprops.csv', index_col=0)
    #
    # features = []
    # for struct in tqdm(concated.loc[:, concated.columns != 'value']['before']):
    #     strl = re.findall('..', struct)
    #     temp = []
    #     for dyad in strl:
    #         temp.append(diprodb[diprodb['PropertyName'] == dyad].values.tolist()[0])
    #     features.append(pd.DataFrame(temp, columns=diprodb.columns).sum())
    # features = pd.DataFrame(features)
    #
    # tsne = TSNE(n_jobs=8)
    # transformed = tsne.fit_transform(features.loc[:, features.columns != 'PropertyName'])
    #
    # vis_x = transformed[:, 0]
    # vis_y = transformed[:, 1]
    #
    # plt.scatter(vis_x, vis_y, c=concated['value'].map({0: 'blue', 1: 'orange'}))
    # # plt.clim(-0.5, 9.5)
    # plt.show()


def plot_feat_import(features, model, y_test, y_pred):
    from sklearn.metrics import roc_auc_score
    df_fi = pd.DataFrame(features.columns[1:], columns=['feature'])
    df_fi['importance'] = list(model.feature_importance('gain'))
    df_fi.sort_values('importance', ascending=False, inplace=True)
    # print(df_fi)
    plt.figure()
    df_fi.head(10).plot(kind='barh', x='feature', y='importance')
    plt.title('Roc_auc is {}'.format(roc_auc_score(y_test, y_pred)))
    plt.xlabel('relative importance')
    return plt


def plot_centered(df, uplifted):
    """
Center uplifted on "peak_pos", create array with number of sec struct in df in range from it. Warning: kind of slow.
    :param df: sec. structs
    :param uplifted: nucleosome positioning
    :return:
    """
    cumsum_arr = np.zeros(1000)
    for peak_pos in tqdm(pd.to_numeric(uplifted.peak_pos)):
        temp_df = df[(df['start'] > peak_pos - 500) & (df['end'] < peak_pos + 500)]
        for ind, row in temp_df.iterrows():
            for i in range(int(row['start']) - peak_pos + 500, int(row['end']) - peak_pos + 500):
                cumsum_arr[i] += 1
    x = np.arange(-500, 500)
    # plt.plot(x, cumsum_arr)
    # plt.ylabel('Количество вторичных структур')
    # plt.xlabel('Позиция относительно центра нуклеосомы')
    # plt.ylim(75000, 110000)
    return x, cumsum_arr


def render_mpl_table(data, title='test', col_width=3.0, row_height=0.625, font_size=14,
                     header_color='#40466e', row_colors=['#f1f1f2', 'w'], edge_color='w',
                     bbox=[0, 0, 1, 1], header_columns=0,
                     ax=None, **kwargs):
    if ax is None:
        size = (np.array(data.shape[::-1]) + np.array([0, 1])) * np.array([col_width, row_height])
        fig, ax = plt.subplots(figsize=size)
        ax.axis('off')
    fig.suptitle(title, fontsize=15)

    mpl_table = ax.table(cellText=data.values, bbox=bbox, colLabels=data.columns, **kwargs)

    mpl_table.auto_set_font_size(False)
    mpl_table.set_fontsize(font_size)
    for k, cell in six.iteritems(mpl_table._cells):
        cell.set_edgecolor(edge_color)
        if k[0] == 0 or k[1] < header_columns:
            cell.set_text_props(weight='bold', color='w')
            cell.set_facecolor(header_color)
        else:
            cell.set_facecolor(row_colors[k[0] % len(row_colors)])
    return ax
