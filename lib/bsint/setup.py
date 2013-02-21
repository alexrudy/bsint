#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 
#  setup.py
#  glbsint
#  
#  Created by Alexander Rudy on 2013-02-12.
#  Copyright 2012 Alexander Rudy. All rights reserved.
# 

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('bsint',parent_package,top_path,
    author = "Alexander Rudy",
    author_email = "arrudy@ucsc.edu",
    )
    config.add_extension('_bsint',['bsint.f','bsint.pyf'],f2py_options=["skip: bsstep pzextr mmid :"])
    print "v",config.get_version()
    return config

if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
    