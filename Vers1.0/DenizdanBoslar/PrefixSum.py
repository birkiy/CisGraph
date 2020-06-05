

import random

array = [random.randint(0, 10**9) for _ in range(3*10**5)]

query = [random.randint(0, 1000) for _ in range(3*10**5)]


L = [0]
tmp = 0
for ar in array:
    
    tmp += ar
    L.append(tmp)



left = 0
for que in query:
    right = que
    s = L[right] - L[left]
    print(s)
