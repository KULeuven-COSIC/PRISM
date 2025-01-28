import json
from utilities.strategy import optimised_strategy_old, optimised_strategy

def write_json(data, file):
    with open(file, 'w', encoding='utf-8') as f: 
        json.dump(data, f, ensure_ascii=False, indent=4)
dim1={122: optimised_strategy_old(122//2, mul_c=2)}
dim2={}
for i in range(1, 256):
    dim2[i]=optimised_strategy(i)
write_json(dim1, 'helpers/strategies_dim1.json')
write_json(dim2, 'helpers/strategies_dim2.json')
