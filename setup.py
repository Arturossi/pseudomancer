from setuptools import find_packages, setup

with open(path.join(this_directory, "requirements.txt"), encoding="utf-8") as f:
    requirements = f.read().splitlines()


def readme():
    with open(path.join(this_directory, "README.md"), encoding="utf-8") as f:
        return f.read()


setup(
    name="pseudomancer",
    version="0.0.1",
    description="pseudomancer: a command line tool for reconstructing "
    + "pseudogenes in prokaryotic genomes",
    long_description=readme(),
    long_description_content_type="text/markdown",
    url="https://github.com/adamd3/pseudomancer",
    license="MIT",
    packages=find_packages(),
    include_package_data=True,
    install_requires=requirements,
    extras_require={"test": "pytest"},
)
