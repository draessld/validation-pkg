from setuptools import setup, find_packages

setup(
    name="validation_pkg",
    version="0.1.0",
    author="Dominika Bohuslavová",
    author_email="dominikadraesslerova@gmail.com",
    description="Simple package to manage input validation for GMO pipeline",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/draessld/validation-pkg",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.7",
    install_requires=[
        # Add your package dependencies here
        # Example: "numpy>=1.19.0",
    ],
)