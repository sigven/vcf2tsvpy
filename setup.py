import setuptools

setuptools.setup(
    name='vcf2tsv',
    version='0.5.0', # versioned by bump2version
    license='MIT',
    author='Sigve Nakken',
    author_email='sigven@ifi.uio.no',
    description='vcf2tsv - conversion of Variant Call Format (VCF) to TSV',
    url='https://github.com/sigven/vcf2tsv',
    packages=setuptools.find_packages(),
    entry_points={
        'console_scripts': 'vcf2tsv = vcf2tsv.main:cli'
    }
)
