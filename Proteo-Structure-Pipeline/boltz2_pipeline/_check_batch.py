import json

data = json.load(open(r'D:\ClawNative0407\vwf-vwd-genomics-research\output\boltz2_blind_scan\batch_001.json'))
print('Batch name:', data['name'])
print('Total jobs:', len(data['jobs']))
print()
for j in data['jobs'][:4]:
    chains = []
    for s in j['sequences']:
        cid = s['protein']['id']
        clen = len(s['protein']['sequence'])
        chains.append(cid + '(' + str(clen) + 'aa)')
    print('  ' + j['name'])
    print('    Chains: ' + ', '.join(chains))
    print('    Properties: ' + str(j['properties']))
    print()
