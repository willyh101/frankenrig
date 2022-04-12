import setuptools

with open('README.md', 'r') as rmf:
    readme = rmf.read()
    
with open('VERSION', 'r') as verf:
    version = verf.read()
    
setuptools.setup(
    name='vispy',
    version=version,
    author='Will Hendricks',
    author_email='hendricksw@berkeley.edu',
    description='Visual stim maker and displayer for 2p imaging in vivo.',
    long_description=readme,
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(),
)