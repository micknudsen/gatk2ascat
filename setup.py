from setuptools import setup, find_packages

setup(

    name='gatk2ascat',
    version='0.1',

    packages=find_packages('src'),
    package_dir={'': 'src'},

    test_suite='tests',

    entry_points={
        'console_scripts': ['gatk2ascat = gatk2ascat.client:main']
    },

    python_requires='>=3.6',

    author='Michael Knudsen',
    author_email='micknudsen@gmail.com',
    license='MIT'

)
