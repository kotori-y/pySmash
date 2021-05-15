# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 13:21:07 2020

@Author: Zhi-Jiang Yang, Dong-Sheng Cao
@Institution: CBDD Group, Xiangya School of Pharmaceutical Science, CSU, China
@Homepage: http://www.scbdd.com
@Mail: yzjkid9@gmail.com; oriental-cds@163.com
@Blog: https://blog.moyule.me

I love my senpai forerver!:P
"""

from __future__ import absolute_import
from distutils.core import setup

# print（__doc__）

package_data= {'smash':['fragments/*']}


setup(name="kotori-smash",  
      version="1.0",
      license="MIT",
      description="",
      long_description="PySmash, a user-friendly and integrated tool which offers convenience and advantage for the automatic structural alerts derivation and application is presented. The current version of PySmash provides both software and Python package, which achieves both flexibility and ease of pipeline integration. Three kinds of substructure derivation algorithm, including circular, path-based and function group-based algorithms are provided, where researchers can customize own parameters about substructure size, accuracy and coverage, statistical significance and server deployment. Besides, PySmash also provides prediction function for derived substructure application.",
      author="Zhi-Jiang Yang (Kotori), Dong-Sheng Cao",
      author_email="yzjkid9@gmail.com",
      maintainer="Zhi-Jiang Yang (Kotori)",
      maintainer_email="kotori@cbdd.me",
      url="https://github.com/kotori-y/pysmash",
      package_data=package_data,
      # include_package_data=True,
      package_dir={'smash':'smash'},
      #py_modules = ['scopy.ScoConfig'],
      packages=['smash']
      )
