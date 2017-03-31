"""
https://archive.ics.uci.edu/ml/datasets/Tic-Tac-Toe+Endgame

https://archive.ics.uci.edu/ml/machine-learning-databases/tic-tac-toe/tic-tac-toe.names

"""
import os

import numpy as np
import tabular as tb

import mine
import utils


din = os.path.join('..', 'data', 'bcancer')
dout = os.path.join('..', 'data', 'CrossValidation')
fdata = os.path.join(din, 'bcancer.data')
ftest = os.path.join(din, 'bcancer.test')
fnames = os.path.join(din, 'bcancer.names')
fcomplete = os.path.join(din, 'bcancer-filtered.csv')
bout = os.path.join('..', 'data', 'bcancer_binary.csv')
fout = os.path.join(din, 'bcancer.csv')

seed = sum([20, 9, 3, 20, 1, 3, 20, 15, 5]) # t:20 i:9 c:3 t:20 a:1 c:3 t:20 o:15 e:5
num_folds = 10
max_cardinality = 3
min_support = 0.05
labels = ['No', 'Yes']
minor = True

np.random.seed(seed)

#x = open(bout, 'rU').read().strip().split('\n')
#x = [line.split(' ') for line in x]
#print x, len(x)

if not os.path.exists(din):
    os.mkdir(din)

if not os.path.exists(dout):
    os.mkdir(dout)

if not os.path.exists(fdata):
    print 'downloading data'
    uroot = 'https://archive.ics.uci.edu/ml/machine-learning-databases/breast-cancer-wisconsin/'
    os.system('wget %sbreast-cancer-wisconsin.data -O %s' % (uroot, fdata))
    os.system('wget %sbreast-cancer-wisconsin.names -O %s' % (uroot, fnames))

names = ['class', 'age', 'menopaus', 'c4', 'c5', 'c6', 'c7', 'c8', 'c9', 'x_wins']

print 'read original train data:', fdata
x = open(fdata, 'rU').read().strip().split('\n')
x = [','.join(line.split(', ')) for line in x if '?' not in line]
assert (len(x) == 958)

f = open(fcomplete, 'w')
f.write(','.join(names) + '\n')
f.write('\n'.join(x))
f.close()

print 'lightly process data (e.g., to make binary features)'
x = tb.tabarray(SVfile=fcomplete)

x_wins = np.cast[int](x['x_wins'] == 'positive')

y = x[names[:-1]].colstack(tb.tabarray(columns=[x_wins], names=names[-1:]))

print 'write categorical dataset', fout
y.saveSV(fout)

print 'write binary dataset', bout
b = utils.to_binary(y)
b.saveSV(bout)


print 'permute and partition dataset'
split_ind = np.split(np.random.permutation(len(x) / num_folds * num_folds), num_folds)
print 'number of folds:', num_folds
print 'train size:', len(split_ind[0]) * (num_folds - 1)
print 'test size:', len(split_ind[0])

num_rules = np.zeros(num_folds, int)
for i in range(num_folds):
    print 'generate cross-validation split', i
    cv_root = 'bcancer_%d' % i
    test_root = '%s_test' % cv_root
    train_root = '%s_train' % cv_root
    ftest = os.path.join(dout, '%s.csv' % test_root)
    ftrain = os.path.join(dout, '%s.csv' % train_root)
    btest = os.path.join(dout, '%s-binary.csv' % test_root)
    btrain = os.path.join(dout, '%s-binary.csv' % train_root)
    train_ind = np.concatenate([split_ind[j] for j in range(num_folds) if (j != i)])
    y[split_ind[i]].saveSV(ftest)
    y[train_ind].saveSV(ftrain)
    b[split_ind[i]].saveSV(btest)
    b[train_ind].saveSV(btrain)

    print 'mine rules from', ftrain
    """
    num_rules[i] = mine.mine_binary(din=dout, froot=train_root,
                                   max_cardinality=max_cardinality,
                                   min_support=min_support,
                                   minor=minor)
    mine.apply_binary(din=dout, froot=cv_root)
    """

    num_rules[i] = mine.mine_rules(din=dout, froot=train_root,
                                   max_cardinality=max_cardinality,
                                   min_support=min_support,
                                   minor=minor, labels=labels)
    mine.apply_rules(din=dout, froot=cv_root, labels=labels)

print '(min, max) # rules mined per fold:', (num_rules.min(), num_rules.max())

#ben.driver(din='../data/tictactoe', dout='../data/tictactoe', froot='tictactoe', train_suffix='.csv',
#           delimiter=',', is_binary=False, maxlhs=2, minsupport=2.5, out_suffix='')
#minority.compute_minority(froot='tictactoe', dir='../data/tictactoe')

#mine.mine_rules(din=din, froot=root, max_cardinality=max_cardinality,
#                min_support=min_support, labels=labels, suffix='_e', minor=minor)

