import pandas as pd

df = pd.read_excel("元素电子结构.xlsx",header=1)
# headers = df.columns.tolist()
# for each in headers:
#     if each in ['idx', 'sym','Unnamed: 61']:
#         continue
#     l,m,n=each.split(",")
#     print(f"{l:>2},{m:>2},{n:>2}")

def get_idx(row:int):
    if row<=2:
        return 0
    if 2<row<=10:
        return 2
    if 10<row<=18:
        return 10
    if 18<row<=36:
        return 18
    if 36<row<=54:
        return 36
    if 54<row<=86:
        return 54
    if 86<row<=118:
        return 86
    assert row <= 118

for idx,row in df.iterrows():
    row_list = row.tolist()
