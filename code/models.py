from sklearn.model_selection import train_test_split
import numpy as np
# %load_ext wurlitzer  # c++ output to jupyter
# import telepyth  # push notif in telegram
# %telepyth - t 1260389131217015787


def lgbm(feats, concated):
    """
fit lgbm with default params
    :param feats: DataFrame to be splitted and fitted
    :param concated: from this Dataframe target is taken
    :return: model and plt for plt.show()
    """
    x_train, x_test, y_train, y_test = train_test_split(
        feats.drop('PropertyName', axis=1).values, concated['value'].values,
        test_size=0.3)


    import lightgbm as lgb
    from sklearn.metrics import roc_auc_score
    np.random.seed(42)
    params = {
    #     "max_bin": 1024,
    #     "learning_rate": 0.01,
        "boosting_type": "goss",
        "objective": "binary",
    #     'num_iterations':10000,
        "metric": "auc",
    #     "num_leaves": 10000,
        "verbose": 1,
    #     "min_data": 100,
    #     "boost_from_average": True,
        'early_stopping_round': 50
    }

    d_train = lgb.Dataset(x_train, y_train)
    d_valid = lgb.Dataset(x_test, label=y_test)
    model = lgb.train(params, d_train, valid_sets=d_valid)
    # model.save_model('lgbm')
    # model = lgb.Booster(model_file='lgbm')

    y_pred = model.predict(x_test)
    print(roc_auc_score(y_test, y_pred))
    from code.visualise import plot_feat_import
    plot = plot_feat_import(feats, model, y_test, y_pred)
    return model, plot