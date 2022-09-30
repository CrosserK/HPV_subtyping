
import difflib

terminal="ReadForCall11_MinLen15_MinQual20_chr17_deletions_12_10_1000MbpData_30000Mean_15000stdev.Sorted.Minimap2_Sniffles_sorted.txt"
adress=""

cases=[(terminal,adress)]

for a,b in cases:
    print('{} => {}'.format(a,b))
    for i,s in enumerate(difflib.ndiff(a, b)):
        if s[0]==' ': continue
        elif s[0]=='-':
            print(u'Delete "{}" from position {}'.format(s[-1],i))
        elif s[0]=='+':
            print(u'Add "{}" to position {}'.format(s[-1],i))
    print()

def sameOrNot(case1,case2):
    if case1 == case2:
        return "Yes"
    else:
        return "No"

print("Are they exactly the same:", sameOrNot(terminal,adress))
