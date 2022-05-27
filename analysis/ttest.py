import pandas, scipy, openpyxl
import scipy.stats

fname = "data/ttest_input.xlsx"

############################################################

df = pandas.read_excel(fname)
colAs = ['SCR1','SCR2']
colB = 'YY1'
res = []
print(f"Value, equal variance, sample 1, smaple 2, P value")
for colA in colAs:
    stat, p = scipy.stats.ttest_ind(df[colA], df[colB])
    print(f"P value,equal variance,{colA},{colB},{p}")
    res.append({"equal variance?": "equal variance",
               "sample 1": colA, 'sample 2': colB,
               'P value': p})
    stat, p = scipy.stats.ttest_ind(df[colA], df[colB], equal_var=False)
    print(f"P value,unequal variance,{colA},{colB},{p}")
    res.append({"equal variance?": "unequal variance",
               "sample 1": colA, 'sample 2': colB,
               'P value': p})
print('-' * 50)

done = []

print(f"Value, equal variance, sample 1, smaple 2, P value")
for colA in colAs:
    for colB in colAs:
        if colA == colB:
            continue
        tup = (colB, colA)
        if tup in done:
            continue
        else:
            done.append((colA, colB))
        stat, p = scipy.stats.ttest_ind(df[colA], df[colB])
        print(f"P value,equal variance,{colA},{colB},{p}")
        res.append({"equal variance?": "equal variance",
                   "sample 1": colA, 'sample 2': colB,
                   'P value': p})
        stat, p = scipy.stats.ttest_ind(df[colA], df[colB], equal_var=False)
        print(f"P value,unequal variance,{colA},{colB},{p}")
        res.append({"equal variance?": "unequal variance",
                   "sample 1": colA, 'sample 2': colB,
                   'P value': p})

res = pandas.DataFrame(res)
print(res)
res.to_excel('results/ttest_output.xlsx')