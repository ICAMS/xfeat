#!/usr/bin/python3

import os
from setuptools import setup, find_packages, Extension

try:
    from Cython.Build import cythonize
except ImportError:
    cythonize = None


def no_cythonize(extensions, **_ignore):
    for extension in extensions:
        sources = []
        for sfile in extension.sources:
            path, ext = os.path.splitext(sfile)
            if ext in (".pyx", ".py"):
                if extension.language == "c++":
                    ext = ".cpp"
                else:
                    ext = ".c"
                sfile = path + ext
            sources.append(sfile)
        extension.sources[:] = sources
    return extensions


extensions = [
    Extension("xfeat.basic", ["src/xfeat/basic.pyx"]),
    Extension("xfeat.main",  ["src/xfeat/main.pyx"]),
]

CYTHONIZE = bool(int(os.getenv("CYTHONIZE", 0))) and cythonize is not None

if CYTHONIZE:
    compiler_directives = {"language_level": 3, "embedsignature": True}
    extensions = cythonize(extensions, compiler_directives=compiler_directives)
else:
    extensions = no_cythonize(extensions)

with open("requirements.txt") as fp:
    install_requires = fp.read().strip().split("\n")

setup(
    ext_modules=extensions,
    install_requires=install_requires,
    extras_require={
        "docs": ["sphinx", "sphinx-rtd-theme"]
    },
)
