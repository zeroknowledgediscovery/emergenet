from setuptools import setup, find_packages
from setuptools.command.install import install
from codecs import open
from os import path
import os, shutil, gzip

package_name = 'emergenet'
example_dir = 'examples/'
bin_dir = 'bin/'
example_data_dir = example_dir + 'example_data/'

version = {}
with open("version.py") as fp:
    exec(fp.read(), version)

here = path.abspath(path.dirname(__file__))

# Get the long description from the relevant file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

class PostInstallCommand(install):
    '''Post-installation for installation mode.'''

    def run(self):
        install.run(self)
        self.unzip_files(path.dirname(__file__))

    def unzip_files(self, directory):
        for root, dirs, files in os.walk(directory):
            for file in files:
                if file.endswith('.gz'):
                    gz_file_path = path.join(root, file)
                    output_file_path = path.join(root, file[:-3])
                    if not path.exists(output_file_path):
                        with gzip.open(gz_file_path, 'rb') as f_in:
                            with open(output_file_path, 'wb') as f_out:
                                shutil.copyfileobj(f_in, f_out)
                        os.remove(gz_file_path)

setup(
    name=package_name,
    author='zed.uchicago.edu',
    author_email='ishanu@uchicago.edu',
    version = str(version['__version__']),
    packages=find_packages(),
    scripts=[],
    url='https://github.com/zeroknowledgediscovery/emergenet',
    license='LICENSE',
    description='Superfast Risk Estimation of Emerging Pathogens',
    keywords=[
        'computational biology',
        'decision trees', 
        'machine learning', 
        'emerging pathogens'],
    download_url='https://github.com/zeroknowledgediscovery/emergenet/archive/'+str(version['__version__'])+'.tar.gz',
    long_description=long_description,
    long_description_content_type='text/x-rst',
    install_requires=[
        "scikit-learn", 
        "scipy", 
        "numpy",  
        "pandas",
        "joblib",
        "quasinet==0.1.54",
        "scipy"],
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 4 - Beta',
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Information Analysis",
        "Topic :: Software Development :: Libraries",
        "Programming Language :: Python :: 3.6"],
    include_package_data=True,
    cmdclass={
        'install': PostInstallCommand,
    },
    )
