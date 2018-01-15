lines = []
with open('install/pip-requirements.txt', 'r') as f:
    lines = f.read() 

with open('install/conda_requirements_anaconda.txt', 'r') as f:
    lines += f.read()
with open('install/conda_requirements_anaconda.txt', 'r') as f:
    lines += f.read()

res = (list('- ' + l  +'\n' for l in lines.split('\n')))

with open ('install/requirements.md', 'wb') as f:
    f.write('# Required Packages beyond miniconda install\n')
    for l in res[:-1]:
        f.write(l)
