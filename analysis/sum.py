import pandas, os, re
df = pandas.read_excel('data/sum_input.xlsx', index_col=None)

df.Protein = [str(x) for x in df.Protein]
df['Rank'] = df.index
df.index = df.Protein
print(df.index)
df.index =  [re.sub('2018-09-02 00:00:00', "'SEPT2", x) for x in df.index]
df.index =  [re.sub('2018-09-07 00:00:00', "'SEPT7", x) for x in df.index]
df.index =  [re.sub('2018-09-11 00:00:00', "'SEPT11", x) for x in df.index]
df.index =  [re.sub('2018-09-09 00:00:00', "'SEPT9", x) for x in df.index]
df.index =  [re.sub('2018-10-04 00:00:00', "'OCT4", x) for x in df.index]
print(df.index)
sub = df.groupby(by='Protein').sum()
sub.to_excel('results/sum_output.xlsx')
