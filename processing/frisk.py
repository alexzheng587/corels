"""
http://www.nyclu.org/content/stop-and-frisk-data


https://5harad.com/papers/frisky.pdf

"we model the likelihood of recovering a weapon in a CPW stop via
logistic regression ... we use only the 18 stop circumstances officers already
consider (listed in Table 1, excluding the two 'other' categories), indicator
variables for each of the 77 precincts and indicator variables for the three
location types (public housing, transit and 'neither'); we do not include
interactions.  To further reduce model complexity and increase interpretability,
we constrain the 18 coefficients corresponding to stop reasons to be non-negative."

Table 1, primary stop circumstance(s):

    Suspicious object, fits description, casing, acting as lookout,
    suspicious clothing, drug transaction, furtive movements, actions
    of violent crime, suspicious bulge and/or other

Table 1, additional stop circumstance(s):

    Witness report, ongoing investigation, proximity to crime
    scene, evasive response, associating with criminals, changed
    direction, high crime area, time of day, sights and sounds of
    criminal activity and/or other

"""
import os

import numpy as np
import tabular as tb

import mine
import utils


city_dict = {1: 'Manhattan', 2: 'Brooklyn', 3: 'Bronx', 4: 'Queens', 5: 'Staten-Island'}

sex_dict = {0: 'female', 1: 'male'}

race_dict = {1: 'black', 2: 'black Hispanic', 3: 'white Hispanic', 4: 'white',
             5: 'Asian/Pacific Islander', 6: 'Am. Indian/Native'}

build_dict = {1: 'heavy', 2: 'muscular', 3: 'medium', 4: 'thin'}

stop_dict = {'cs_objcs': 'reason for stop - suspicious object',
             'cs_descr': 'reason for stop - fits description',
             'cs_casng': 'reason for stop - casing',
             'cs_lkout': 'reason for stop - acting as lookout',
             'cs_cloth': 'reason for stop - suspicious clothing',
             'cs_drgtr': 'reason for stop - drug transaction',
             'cs_furtv': 'reason for stop - furtive movements',
             'cs_vcrim': 'reason for stop - actions of violent crime',
             'cs_bulge': 'reason for stop - suspicious bulge',
             'cs_other': 'reason for stop - other',
             'ac_proxm': 'additional circumstances - proximity to crime scene',
             'ac_evasv': 'additional circumstances - evasive response',
             'ac_assoc': 'additional circumstances - associating with criminals',
             'ac_cgdir': 'additional circumstances - changed direction',
             'ac_incid': 'additional circumstances - high crime area',
             'ac_time' : 'additional circumstances - time of day',
             'ac_stsnd': 'additional circumstances - sights and sounds of criminal activity',
             'ac_rept' : 'additional circumstances - witness report',
             'ac_inves': 'additional circumstances - ongoing investigation',
             'ac_other': 'additional circumstances - other'}

inout_dict = {0: 'outside', 1: 'inside'}
trhsloc_dict = {0: 'neither-housing-nor-transit-authority', 1: 'housing-authority', 2: 'transit-authority'}

def rename_pos(s):
    return s.replace(' - ', '=').replace('reason for stop', 'stop-reason').replace('additional circumstances', 'circumstances').replace(' ', '-')

def rename_neg(s):
    return s.replace(' - ', '=not-').replace('reason for stop', 'stop-reason').replace('additional circumstances', 'circumstances').replace(' ', '-')

def age_func(a):
    if (a < 18):
        return '<18'    # support = 5630
    elif (a <= 21):
        return '18-21'  # support = 9960
    elif (a <= 25):
        return '22-25'  # support = 8044
    elif (a <= 30):
        return '26-30'  # support = 6222
    elif (a <= 40):
        return '31-40'  # support = 7115
    else:
        return '>40'    # support = 6776


din = os.path.join('..', 'data', 'frisk')
dout = os.path.join('..', 'data', 'CrossValidation')
zdata = os.path.join('..', 'data', '2014-20SQF.zip')
fdata = os.path.join(din, '2014-SQF-web.csv')
fout = os.path.join(din, 'frisk.csv')
bout = os.path.join('..', 'data', 'frisk-binary.csv')
fout = os.path.join(din, 'frisk.csv')

seed = sum([1, 4, 21, 12, 20]) # f:6, r:18, i:09, s:19, k:11
num_folds = 10
max_cardinality = 1
min_support = 0.001
labels = ['no', 'yes']
minor = True
predict_frisked = False
predict_weapon = True

np.random.seed(seed)

if not os.path.exists(din):
    os.mkdir(din)

if not os.path.exists(dout):
    os.mkdir(dout)

if not os.path.exists(fdata):
    print 'downloading data'
    uroot = 'http://www.nyclu.org/files/2014%20SQF.zip'
    os.system('wget %s -O %s' % (uroot, zdata))
    print 'unzipping data'
    os.system('unzip %s' % zdata)
    os.system('mv 2014\ SQF/* %s' % din)
    print 'renaming files'
    for f in os.listdir(din):
        f1 = os.path.join(din, f)
        f2 = f1.replace(' ', '-').replace('_', '-')
        f1 = f1.replace(' ', '\\ ')
        os.system('mv %s %s' % (f1, f2))

print 'lightly process data (e.g., to make binary features)'
x = tb.tabarray(SVfile=fdata)
assert (len(x) == 45787)

assert not np.isnan(x['frisked']).any()
assert x['frisked'].sum() == 30345  # 66% frisked

assert not np.isnan(x['searched']).any()
assert x['searched'].sum() == 7283  # 16% searched

assert (x['frisked'] & x['searched']).sum() == 6667
assert (x['searched'] & np.invert(x['frisked'])).sum() == 616

weapon_list = ['pistol', 'riflshot', 'asltweap', 'knifcuti', 'machgun', 'othrweap']
w = x[weapon_list].extract()
(ii, jj) = np.isnan(w).nonzero()
w[ii, jj] = 0
weapon = w.any(axis=1)
assert weapon.sum() == 1520
assert (x['searched'] & weapon).sum() == 1081   # 15% of searches yield weapon
assert (x['frisked'] & weapon).sum() == 1414    # 4.7% of frisks lead to weapon (possibly via search)
assert (x['searched'] | x['frisked']).sum() == 30961
assert ((x['searched'] | x['frisked']) & weapon).sum() == 1445

assert not np.isnan(x['arstmade']).any()
assert x['arstmade'].sum() == 6898  # 15% arrested

assert len(set(x['year'])) == 1     # year of stop

assert len(set(x['pct'])) == 77     # precinct of stop
assert (x['pct'].min() >= 1) and (x['pct'].max() <= 123)

assert len(set(x['ser_num'])) == 2281   # UF-250 serial number

# datestop, timestop

assert set(x['city']) == set(range(1, 6))

assert np.isnan(x['sex']).sum() == 394
assert np.isnan(x['race']).sum() == 1039

# dob

assert np.isnan(x['age']).sum() == 107
assert (x['age'] >= 100).sum() == 27
assert x[x['age']==366]['dob']==6161978
x[x['age']==366]['age'] = 36
assert (x['age'] >= 90).sum() == 29     # the two entries in the 90s are ambiguous (see dob)
assert (x['age'] < 12).sum() == 99      # the single digit / younger ages seem like typos (see height, weight)

assert not np.isnan(x['height']).any()
assert not np.isnan(x['weight']).any()

assert np.isnan(x['haircolr']).sum() == 492
assert np.isnan(x['eyecolor']).sum() == 286
assert np.isnan(x['build']).sum() == 721

assert len(set(x['othfeatr'])) == 220

assert set(x['inout']) == set([0, 1])
assert set(x['trhsloc']) == set([0, 1, 2])

stop_reasons = [n for n in x.dtype.names if (n.startswith('cs') or n.startswith('ac'))
                                      and ('other' not in n)] # or n.startswith('ac'))
for sr in stop_reasons:
    x[sr][np.isnan(x[sr]).nonzero()[0]] = 0

names = ['city', 'sex', 'race', 'age', 'build']

keep = np.invert(np.isnan(x[names].extract()).any(axis=1))
x = x[keep]
weapon = weapon[keep]
assert len(x) == 43858  # throw out 1929 records with missing data
ikeep = (x['age'] > 11) & (x['age'] < 90)
x = x[ikeep]     # throw out age extremes
weapon = weapon[ikeep]
assert len(x) == 43747

if (predict_weapon):
    ikeep = (x['searched'] == 1) | (x['frisked'] == 1)
    x = x[ikeep]
    weapon = np.cast[int](weapon[ikeep])
    assert (len(x) == 29595)

city = [city_dict[i] for i in x['city']]
sex = [sex_dict[i] for i in x['sex']]
race = [race_dict[i] for i in x['race']]
age = [age_func(i) for i in x['age']]
build = [build_dict[i] for i in x['build']]
inout = [inout_dict[i] for i in x['inout']]
location = [trhsloc_dict[i] for i in x['trhsloc']]

stop_reasons_list = []
for sr in stop_reasons:
    stop_reasons_list += [[rename_pos(stop_dict[sr]) if s else rename_neg(stop_dict[sr]) for s in x[sr]]]

if (predict_frisked):
    columns = [city, sex, race, age, build, inout, location] + stop_reasons_list + [x['frisked']]
    cnames = names + ['inout', 'location'] + stop_reasons + ['frisked']
elif (predict_weapon):
    columns = stop_reasons_list + [city, location, inout, weapon]
    cnames = stop_reasons + ['city', 'location', 'inout', 'weapon']

print 'write categorical dataset', fout
y = tb.tabarray(columns=columns, names=cnames)

if (predict_weapon):
    y0 = y[y['weapon'] == 0]
    y1 = y[y['weapon'] == 1]
    y = y0.rowstack(y1)

y.saveSV(fout)

print 'write binary dataset', bout
b = utils.to_binary(y)
b.saveSV(bout)

if (not predict_weapon):
    split_ind = np.split(np.random.permutation(len(y) / num_folds * num_folds), num_folds)
else:
    s0 = np.split(np.random.permutation(len(y0) / num_folds * num_folds), num_folds)
    s1 = np.split(len(y0) + np.random.permutation(len(y1) / num_folds * num_folds), num_folds)
    s1 = [i1[np.random.randint(0, len(i1), len(s0[0]))] for i1 in s1]
    split_ind = [np.concatenate([i0, i1]) for (i0, i1) in zip(s0, s1)]

print 'permute and partition dataset'
print 'number of folds:', num_folds
print 'train size:', len(split_ind[0]) * (num_folds - 1)
print 'test size:', len(split_ind[0])

num_rules = np.zeros(num_folds, int)
for i in range(num_folds):
    print 'generate cross-validation split', i
    cv_root = 'frisk_%d' % i
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
    num_rules[i] = mine.mine_rules(din=dout, froot=train_root,
                                    max_cardinality=max_cardinality,
                                    min_support=min_support, labels=labels,
                                    minor=minor)
    mine.apply_rules(din=dout, froot=cv_root, labels=labels)

print '(min, max) # rules mined per fold:', (num_rules.min(), num_rules.max())
