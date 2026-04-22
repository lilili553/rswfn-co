import os
from pathlib import Path
from openai import OpenAI

text=Path("26.04.01/DFTxclib.F").read_text()
lines=text.split("\n")
marks=[]
names=[]
for l,line in enumerate(lines):
    if "subroutine" in line:
        marks.append(l)
        loc=line.index("subroutine")+11
        names.append(line[loc:])
marks.append(len(lines))

texts=[]
for i in range(len(marks)-1):
    up=marks[i]
    lo=marks[i+1]
    text="\n".join(lines[up:lo])
    texts.append(text)
print(names)
# print(texts[0])
for name,text in zip(names,texts):
    print(name)
    Path(f"26.04.01/{name}.f90").write_text(text)