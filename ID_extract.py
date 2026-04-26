lines=open('filtered_nucl.txt').readlines()
ids=[l.split('\t')[0].strip() for l in lines]
open('ids.txt','w').write('\n'.join(ids))
print(len(ids),'IDs saved')
