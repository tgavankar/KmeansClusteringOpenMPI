import csv
import re
import os
import glob
 
path = 'large_results/cProfile/'
writer = csv.writer(open("analysis.csv", "w"))

for infile in glob.glob( os.path.join(path, '*.txt') ):
    txt=infile
    re1='((?:[a-z][a-z]+))'	# Word 1
    re2='(\\d+)'	# Integer Number 1
    re3='.*?'	# Non-greedy match on filler
    re4='(\\d+)'	# Integer Number 2
    re5='.*?'	# Non-greedy match on filler
    re6='(\\d+)'	# Integer Number 3
    re7='((?:[a-z][a-z]+))'	# Word 2

    rg = re.compile(re1+re2+re3+re4+re5+re6+re7,re.IGNORECASE|re.DOTALL)
    m = rg.search(txt)
    if m:
        o = []
        o.append(m.group(1).strip())
        o.append(m.group(2).strip())
        o.append(m.group(3).strip())
        o.append(m.group(4).strip())
        o.append(m.group(5).strip())
       
    
    with open(infile, 'r') as f:
        re1='.*?'	# Non-greedy match on filler
        re2='([+-]?\\d*\\.\\d+)(?![-+0-9\\.])'	# Float 1

        rg = re.compile(re1+re2,re.IGNORECASE|re.DOTALL)
        of = []
        for txt in f:
            if txt[-8:-1] != "seconds":
                continue
            m = rg.search(txt)
            if m:
                of.append(m.group(1))
                
        o.append(of[-1].strip())

    writer.writerow(o)
