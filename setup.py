#!/usr/bin/env python
# -*- coding: utf-8 -*-

# {# pkglts, pysetup.kwds
# format setup arguments

from setuptools import setup, find_packages


short_descr = "Python library to run biological simulations of 3D, multicellular tissues."
readme = open('README.rst').read()
history = open('HISTORY.rst').read().replace('.. :changelog:', '')


# find version number in src/multicell/version.py
version = {}
with open("src/multicell/version.py") as fp:
    exec(fp.read(), version)


setup_kwds = dict(
    name='multicell',
    version=version["__version__"],
    description=short_descr,
    long_description=readme + '\n\n' + history,
    author="Jean-Louis Dinh, ",
    author_email="contact@jldinh.com, ",
    url='https://github.com/Jean-Louis Dinh/multicell',
    license='mit',
    zip_safe=False,

    packages=find_packages('src'),
    package_dir={'': 'src'},
    install_requires=[
        ],
    tests_require=[
        "coverage",
        "mock",
        "nose",
        "twine",
        ],
    entry_points={},
    keywords='',
    
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python",
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 2 :: Only",
        "Programming Language :: Python :: 2.7",
    ],
    test_suite='nose.collector',
)

# #}

# change setup_kwds below before the next pkglts tag

# do not change things below
# {# pkglts, pysetup.call
setup(**setup_kwds)

# #}

