# %%
import os
import glob

# %%
fqfile = list(map(lambda x: os.path.basename(x), glob.glob('**/*.fastq.gz', recursive=True)))
print('samples:')
samples = set(map(lambda x: x.split('_')[0], fqfile))
for z in samples:
    print(' '*4, f'{z}:')


