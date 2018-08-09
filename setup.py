from setuptools import setup
from numpy.distutils.core import setup, Extension

VERSION = '0.1.0'

tessmaps = Extension("src._tessmaps",
                    ["src/tessmaps.py",
                     "src/get_time_on_silicon.py",
                     "src/plot_skymaps_of_targets.py",
                     "src/get_data.py",
                     "src/get_targets.py"
                     ])

setup(
    name='tessmaps',
    version=VERSION,
    description="some tools for understanding what TESS is looking at",
    author='Luke Bouma',
    author_email='luke@astro.princeton.edu',
    url='https://github.com/lgbouma/tessmaps',
    install_requires=['numpy', 'astropy', 'astroquery', 'matplotlib', 'scipy',
                      'seaborn'],
    packages=['tessmaps', 'tessmaps.data'],
    package_dir={'tessmaps': 'src', 'tessmaps.data': 'data'},
    package_data={'tessmaps.data': ['*.csv', '*.vot']},
    ext_modules=[tessmaps,],
    classifiers=[
          "Development Status :: 3 - Alpha",
          "License :: OSI Approved :: MIT License",
          "Operating System :: OS Independent",
          "Programming Language :: Python",
          "Intended Audience :: Science/Research",
          "Topic :: Scientific/Engineering :: Astronomy", ],
    include_package_data=True
)
