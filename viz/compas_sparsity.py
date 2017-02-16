import matplotlib.patches as mp
import matplotlib.pyplot as plt
import numpy as np
import tabular as tb

with_training = False

plt.ion()

z = tb.tabarray(SVfile='../eval/compas_sparsity_train.csv')
y = tb.tabarray(SVfile='../eval/compas_sparsity-CORELS.csv', names=z.dtype.names[:-1], namesinheader=False)

x = z[z['Method'] != 'CORELS'].rowstack(y[y['Method'] == 'CORELS'])

m = x.aggregate(On=['Method', 'C', 'cp', 'R'], AggFuncDict={'accuracy': np.mean, 'leaves': np.mean, 'train_accuracy': np.mean})
s = x.aggregate(On=['Method', 'C', 'cp', 'R'], AggFuncDict={'accuracy': np.std, 'leaves': np.std, 'train_accuracy': np.std})

fig = plt.figure(1, figsize=(8, 3.5))
plt.clf()
ax = plt.subplot2grid((20, 1), (0, 1), colspan=1, rowspan=19)

m.sort(order=['Method', 'C', 'cp', 'R'])
s.sort(order=['Method', 'C', 'cp', 'R'])

ind = range(10, 15) + range(5, 10) + range(5)
m = m[ind].copy()
s = s[ind].copy()

#m = m[m['Method'] == 'CORELS'].rowstack(m[m['Method'] != 'CORELS'])
#.rowstack(m[m['Method'] == 'RIPPER']).rowstack(m[m['Method'] == 'SBRL']).rowstack(m[m['Method'] == 'CART']).rowstack(m[m['Method'] == 'C4.5'])
#s = s[s['Method'] == 'CORELS'].rowstack(s[s['Method'] != 'CORELS'])
#.rowstack(s[s['Method'] == 'RIPPER']).rowstack(s[s['Method'] == 'SBRL']).rowstack(s[s['Method'] == 'CART']).rowstack(s[s['Method'] == 'C4.5'])

data = zip(m['Method'], m['leaves'], m['accuracy'],  s['leaves'], s['accuracy'], m['train_accuracy'], s['train_accuracy'])

ms = 5
cdict = {'CORELS': 'r', 'C4.5': 'c', 'CART': 'gray', 'RIPPER': 'k', 'SBRL': 'purple'}
mdict = {'CORELS': 's', 'C4.5': '^', 'CART': 'd', 'RIPPER': 'v', 'SBRL': 'o'}
msdict = {'CORELS': 10, 'C4.5': ms, 'CART': ms, 'RIPPER': ms*2, 'SBRL': ms*2}
mfcdict = {'CORELS': 'coral', 'C4.5': 'paleturquoise', 'CART': 'white', 'RIPPER': 'lightgray', 'SBRL': 'plum'}
msvec = np.array([12, 8, 4, 10, 10, 12, 10, 8, 6, 4, 4, 6, 8, 10, 12])
mew = 2

i = 0
for (method, xx, yy, w, h, ty, th) in data:
    if (0):#(i == 8):
        mfc = 'None'
    else:
        mfc = mfcdict[method]
    if (w == 0):
        plt.plot(xx, yy, color=cdict[method], linewidth=0, marker=mdict[method], markersize=msvec[i], markeredgewidth=mew, markeredgecolor=cdict[method], markerfacecolor=mfc)
    else:
        plt.errorbar(xx, yy, xerr=w, color=cdict[method], linewidth=0, marker=mdict[method], markersize=msvec[i], markeredgewidth=mew, markeredgecolor=cdict[method], markerfacecolor=mfc, capsize=4, elinewidth=2)
    i += 1

i = 0
for (method, xx, yy, w, h, ty, th) in data:
    if (0):#(i == 8):
        mfc = 'None'
    else:
        mfc = mfcdict[method]
    if (i == 2):
        xx -= 0.075
    plt.errorbar(xx, yy, yerr=h, color=cdict[method], linewidth=0, marker=mdict[method], markersize=msvec[i], markeredgewidth=mew, markeredgecolor=cdict[method], markerfacecolor=mfc, capsize=0, elinewidth=1)
    i += 1


if (with_training):
    i = 0
    for (method, xx, yy, w, h, ty, th) in data:
        if (i == 8):
            mfc = 'None'
        else:
            mfc = mfcdict[method]
        if ty:
            plt.plot(xx, ty, '.', color=cdict[method])
            plt.plot([xx, xx], [ty, yy], color=cdict[method])
        i += 1

legend = []
for r in m:
    descr = r['Method']
    if r['C']:
        descr += ' (%s)' % ('%1.3f' % r['C']).strip('0')
    elif r['cp']:
        descr += ' (%s)' % ('%1.3f' % r['cp']).strip('0')
    elif r['R']:
        descr += ' (%s)' % ('%1.3f' % r['R']).strip('0')
    legend += [descr]

fs = 14
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
plt.xlabel('Model size', fontsize=fs)
plt.ylabel('Accuracy', fontsize=fs)
plt.legend(legend, loc='lower right', fontsize=fs-3, numpoints=1, ncol=3)

ax.set_xlim(0, 36)

if (with_training):
    ax.set_ylim(0.635, 0.685)
    plt.show()
    plt.savefig('../figs/compas-sparsity-training.pdf')
else:
    ax.set_ylim(0.60, 0.7)
    plt.show()
    plt.savefig('../figs/compas-sparsity.pdf')
