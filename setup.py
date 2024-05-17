from setuptools import setup, find_packages



setup(
    name='BioML',
    version='0.1.0',
    description='A package for DNA, RNA and protein sequence manipulation',
    author='Neermita Bhattacharya',
    author_email='nemowbio@gmail.com',
    url='https://github.com/Neermita18/BioML',  # Replace with your package URL
    packages=find_packages(),
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
   
        
    
    keywords='DNA RNA protein bioinformatics sequence computational biology genetics genome chromosomes inheritance heredity',
    python_requires='>=3.6',
    install_requires=[],
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
)