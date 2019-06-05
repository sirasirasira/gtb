import sys
import numpy as np
import argparse

# parameters
parser = argparse.ArgumentParser(
		prog = 'gradient_boosting',
		description = 'description',
		epilog = 'end',
		add_help = True
		)

parser.add_argument('-t', help ='num of trees', type = int, default = 300)
parser.add_argument('-s', help ='shrinkage', type = float, default = 1.0)
parser.add_argument('-n', help ='needed impurity decrease', type = float, default = 0.0)
parser.add_argument('-d', help ='max depth', type = int, default = 1000000)
parser.add_argument('training_data_file', help ='training data file')
parser.add_argument('test_data_file', help ='test data file')
parser.add_argument('feature_size', help ='feature size', type = int)

args = parser.parse_args()

#print(args)


# make feature vector
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

X_train, y_train = get_data(args.training_data_file, args.feature_size)
X_test, y_test = get_data(args.test_data_file, args.feature_size)


# learning
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score

# model = GradientBoostingClassifier
'''
model = GradientBoostingClassifier(learning_rate = args.s, n_estimators = args.t, max_depth = args.d, min_impurity_decrease = args.n)
model.fit(X_train, y_train)
y_pred = model.predict(X_test)
print("ACC (gradient boosting):", accuracy_score(y_pred, y_test))
'''

# model = RandomForestClassifier
model = RandomForestClassifier(n_estimators = 100, max_depth = 5)
model.fit(X_train, y_train)
y_pred = model.predict(X_test)
print("ACC (random forest):", accuracy_score(y_pred, y_test))

