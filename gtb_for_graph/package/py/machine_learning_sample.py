import sys
import numpy as np

def get_data(fname, m1):
    with open(fname, "r") as f:
        lns = f.readlines()
    m0 = -1
    out = []
    for line in lns:
        token = line.strip().split(':')
        no, y = [int(i) for i in token[0].split()]
        if no > m0: 
            m0 = no
        idx = [int(i) for i in token[1].split()]
        out.append([no, y, idx])
    ret = np.zeros((m0, m1))
    yret = np.zeros(m0)
    for val in out:
        i, y, arr = val
        yret[i-1] = y
        for j in arr:
            ret[i-1, j-1] = 1
    print("data size:", ret.shape)
    return ret, yret

fname_tr = sys.argv[1]
fname_te = sys.argv[2]
feat_siz = sys.argv[3]

X_train, y_train = get_data(fname_tr, int(feat_siz))
X_test, y_test = get_data(fname_te, int(feat_siz))

from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
from sklearn.model_selection import cross_val_score

#X_train, X_test, y_train, y_test = train_test_split(xx, yy, test_size=0.33, random_state=0)

from sklearn.linear_model import LogisticRegression

model = LogisticRegression()
model.fit(X_train, y_train)
y_pred = model.predict(X_test)

print("ACC (logreg): ", accuracy_score(y_pred, y_test))

# res = cross_val_score(model, xx, yy, cv=4)
# print("ACC (logreg): 4-fold cv", np.mean(res), np.std(res))

from sklearn.ensemble import RandomForestClassifier

model = RandomForestClassifier(n_estimators=50)
model.fit(X_train, y_train)
y_pred = model.predict(X_test)

print("ACC (randomforest):", accuracy_score(y_pred, y_test))

# res = cross_val_score(model, xx, yy, cv=4)
# print("ACC (randomforest): 4-fold cv", np.mean(res), np.std(res))

