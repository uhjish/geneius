from setuptools import setup, find_packages, Extension
module1 = Extension('libgeneius.seq._twobit',
                    sources = ['libgeneius/seq/_twobit.c'])


setup(

    name='geneius',
    version='0.0.1',
    package_dir={'libgeneius':'libgeneius','libgeneius.seq':'libgeneius/seq'},
    packages=['libgeneius','libgeneius.seq'],
    py_modules=['simple_geneius'],
    scripts=['geneius_cl.py'],
    ext_modules = [module1]
    )

