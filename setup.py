# coding utf8
import setuptools
from yxquantgene.versions import get_versions

with open('README.md') as f:
    LONG_DESCRIPTION = f.read()

setuptools.setup(
    name="yxquantgene",
    version=get_versions(),
    author="Yuxing Xu",
    author_email="xuyuxing@mail.kib.ac.cn",
    description="Xu Yuxing's personal quantitative genomic tools",
    long_description=LONG_DESCRIPTION,
    long_description_content_type='text/markdown',
    url="https://github.com/SouthernCD/yxquantgene",
    include_package_data = True,

    # entry_points={
    #     "console_scripts": ["HugeP2G = hugep2g.cli:main"]
    # },    

    packages=setuptools.find_packages(),

    install_requires=[
        'scikit-learn',
        'pandas',
        'yxutil',
        'yxseq',
        'statsmodels',
        'matplotlib',
        'pysam',
        'yxmath',
        'numpy',
        'joblib',
        'seaborn',
        'scipy',
        'cyvcf2',
    ],

    python_requires='>=3.5',
)