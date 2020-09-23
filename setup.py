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


setup(name="smash",  
      version="1.0",
      license="MIT",
      description="",
      long_description="",
      author="Zhi-Jiang Yang (Kotori), Dong-Sheng Cao",
      author_email="yzjkid9@gmail.com",
      maintainer="Zhi-Jiang Yang (Kotori)",
      maintainer_email="kotori@cbdd.me",
      url="",
      package_data=package_data,
      # include_package_data=True,
      package_dir={'smash':'smash'},
      #py_modules = ['scopy.ScoConfig'],
      packages=['smash']
      )
