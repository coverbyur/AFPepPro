from distutils.core import setup

with open('README.md') as f:
    readme = ''.join(f.readlines())

#with open('LICENSE') as f:
#    license = f.read()

PACKAGES = ['AFPeptidePro']

setup(
    name='AFPeptidePro',
    version='0.0.1',
    description='Anti-fouling peptide design and selection algorithm',
    long_description=readme,
    author='Clyde Overby',
    author_email='coverby@ur.rochester.edu',
    packages= PACKAGES,
    entry_points=
    {
        'console_scripts': ["AFPepPro=AFPeptidePro.AFPepPro:start"]
    },
    url='https://github.com/coverby/AFPeptidePro',
#    license=license,
)