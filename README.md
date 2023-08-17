# CSA-MEM
## Circular FM-Index


## Most Significant Common Subsequence Chain of MEMs
Pseudo-Code:
```python
L = NULL
i = 0
for sequence in sequences:
    for prev_MEM, MEM in zip(sequence, sequence[1:]):
        b1 = {prev_MEM.query, prev_MEM.query + prev_MEM.length}
        b2 = {MEM.query, MEM.query + MEM.length}
        if intersection(b1,b2):
            remove_smaller(prev_MEM, MEM)
    sort(sequence, reference)
    for prev_MEM, MEM in zip(sequence, sequence[1:]):
        b1 = {prev_MEM.reference, prev_MEM.reference + prev_MEM.length}
        b2 = {MEM.reference, MEM.reference + MEM.length}
        if intersection(b1,b2):
            remove_smaller(prev_MEM, MEM)
    if i==0:
        L = S[0]
    else
        L = L + sequence
        Temp_L = NULL
        for prev_MEM, MEM in zip(L, L[1:]):
            b1 = {prev_MEM.reference, prev_MEM.reference + prev_MEM.length}
            b2 = {MEM.reference, MEM.reference + MEM.length}
            if intersection(b1,b2):
                Temp_L = Temp_L + intersection(b1,b2)
        L = Temp_L
                
    i = i + 1

return longest_chain(L)
```

