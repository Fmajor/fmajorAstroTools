from __future__ import (absolute_import, division, print_function, unicode_literals)
from setuptools import setup, find_packages
from distutils.command.install import install
from distutils.command.build import build
from subprocess import call
import importlib
import os
import sys

# setups
name="fmajorAstroTools"
version="0.0"
packages = ["fmajorAstroTools"]
install_requires=['ipython', 'matplotlib', 'numpy', 'astropy', "colorama", "ipdb", "scipy", "pyds9"]
author='Jin Wu'
author_email='wujinnnnn@.gmail.com'
license='MIT'
description="Astro Tools."
long_description=""

# for scripts
binExcludeExt = ["swp", "pyc"]
binExcluldeFile = []
binToDelete = ["__pycache__"]
binFolders = ["bin", "fmajorAstroTools/DS10/bin"]
scripts = []
for eachBinFolder in binFolders:
    binFolder = os.path.join(os.path.dirname(__file__), eachBinFolder)
    aux = os.listdir(binFolder)
    if len(binToDelete):
        for eachToDelete in binToDelete:
            todoPath = os.path.join(eachBinFolder, eachToDelete)
            if os.path.exists(todoPath):
                os.system("rm -rf {}".format(todoPath))
    scripts.extend([os.path.join(eachBinFolder, each) for each in aux if each[0]!="." and each.split(".")[-1] not in binExcludeExt and each not in binExcluldeFile])
print("scripts to install:\n\t{}".format(scripts))

#<=== for makefile
## in each packages' __init__.py
#__makefile__ = {"./clib/makefile":[]}
#__installFiles__ = ["clib/*.o"]
package_data={}
makefiles=[]
for eachPackage in packages:
    thisModule = importlib.import_module(eachPackage)
    packageRoot = eachPackage
    if hasattr(thisModule, "__makefile__"):
        thisMakeFiles = []
        assert isinstance(thisModule.__makefile__, dict), type(thisModule.__makefile__)
        for i,eachMakefileName in enumerate(thisModule.__makefile__):
            thisMakeFiles.append({
                "makefileName":os.path.join(packageRoot, eachMakefileName),
                "args":thisModule.__makefile__[eachMakefileName]})
        makefiles.extend(thisMakeFiles)
    if hasattr(thisModule, "__installFiles__"):
        package_data[str(packageRoot)] = thisModule.__installFiles__

#print(package_data)
def compile():
    for eachMakefile in makefiles:
        todo = ["make", "-f" ]
        todo.append(os.path.basename(eachMakefile["makefileName"]))
        todo.extend(eachMakefile["args"])
        todo.extend(['CFLAGS=-std=c99 -lm'])
        print("make {}".format(eachMakefile["makefileName"]))
        res = call(todo, cwd=os.path.dirname(eachMakefile["makefileName"]))
        if res!=0:
            print("make error!!")
            sys.exit(1)

class costumBuild(build):
    def run(self):
        self.execute(compile, [], "Compiling extinsions")
        build.run(self)

class costumInstall(install):
    def run(self):
        install.run(self)

def readme():
    with open('README.rst') as f:
        return f.read()

#==> for make file

#python setup.py install --record files.txt
#cat files.txt | xargs rm -rf

setup(name=name,
      version=version,
      description=description,
      long_description=long_description,
      install_requires=install_requires,
      author=author,
      author_email=author_email,
      license=license,
      packages=packages,
      include_package_data=True,
      zip_safe=False,
      cmdclass={'build': costumBuild, "install": costumInstall},
      package_data=package_data,
      scripts=scripts,
      #url='http://do-not-exist.com',
      #keywords='funniest joke comedy flying circus',
      #entry_points = { 'console_scripts': ['binName=Module.file:function'], }
      #dependency_links=['http://github.com/user/repo/tarball/master#egg=package-1.0'],
      #classifiers=['see: http://pypi.python.org/pypi?%3Aaction=list_classifiers'],
)
