import setuptools

setuptools.setup(
    name='vcf2tsvpy',
    version='0.5.0', # versioned by bump2version
    license='MIT',
    author='Sigve Nakken',
    author_email='sigven@ifi.uio.no',
    description='vcf2tsvpy - conversion of Variant Call Format (VCF) to TSV',
    url='https://github.com/sigven/vcf2tsvpy',
    packages=setuptools.find_packages(),
    entry_points={
        'console_scripts': 'vcf2tsvpy = vcf2tsvpy.main:cli'
    }
)
