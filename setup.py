from setuptools import setup, find_packages

setup(
    name="gist",
    version="0.2.0",
    packages=find_packages(),
    install_requires=["dacite"],
    entry_points={
        "console_scripts": [
            "gist = gist.scripts.main:cli",
        ],
    },
)
