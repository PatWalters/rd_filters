from setuptools import setup, find_packages


setup(
    name="rd_filters",
    version="0.1",
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'rd_filters=rd_filters.rd_filters:main',
        ],
    },
    install_requires=['pandas', 'docopt'],
    include_package_data=True,
)
