from setuptools import find_packages, setup
import io
import os
import re

with open(
    os.path.join(this_directory, "requirements.txt"), encoding="utf-8"
) as f:
    requirements = f.read().splitlines()


def readme():
    with open(os.path.join(this_directory, "README.md"), encoding="utf-8") as f:
        return f.read()


def read(*names, **kwargs):
    with io.open(
        os.path.join(os.path.dirname(__file__), *names),
        encoding=kwargs.get("encoding", "utf8"),
    ) as fp:
        return fp.read()


def find_version(*file_paths):
    version_file = read(*file_paths)
    version_match = re.search(
        r"^__version__ = ['\"]([^'\"]*)['\"]", version_file, re.M
    )
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")


setup(
    name="pseudomancer",
    version=find_version("pseudomancer/__init__.py"),
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
