# CSA-MEM
## Circular FM-Index

<img width="1297" alt="image" src="https://github.com/andre99salgado/CSA-MEM/assets/61476040/2b3db5b8-3e75-405b-b810-ee5d8d04978a">
<img width="1323" alt="image" src="https://github.com/andre99salgado/CSA-MEM/assets/61476040/a2fd7f95-ef86-4ccf-abc1-413ed4eccf93">
<img width="1302" alt="image" src="https://github.com/andre99salgado/CSA-MEM/assets/61476040/a5be80e0-ba15-4c09-bc26-09bcec12acd4">
<img width="1300" alt="image" src="https://github.com/andre99salgado/CSA-MEM/assets/61476040/8c78ea56-0841-4876-93b5-b940476b3c18">

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
        sort(L, reference)
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

