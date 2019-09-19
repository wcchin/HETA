from setuptools import setup, find_packages

setup(
    name="heta",

    version="1.0.0",

    author="Benny Chin, Huang C. Y.",
    author_email="wcchin.88@gmail.com",

    packages=['heta'],

    include_package_data=True,

    #url="https://github.com/wcchin/pyreveal",

    license="LICENSE",
    description="HETA: Hierarchical Edge Type Analysis",

    long_description=open("README.md").read(),

    classifiers=[
        'Development Status :: 4 - Beta',

        'Intended Audience :: Developers',
        'Intended Audience :: Information Technology',
        'Intended Audience :: Education',
        'Topic :: Documentation',

         'License :: OSI Approved :: MIT License',

        'Programming Language :: Python :: 3.6',
    ],

    keywords='graph theory, global bridge, local bridge, hierarchical, bond, silk, strength of edge, graph, complex network',

    install_requires=[
        "networkx",
        "numpy",
        "scipy",
        "matplotlib",
        "tqdm",
        "pathos"
    ],
)
