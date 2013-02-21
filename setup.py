#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 
#  setup.py
#  glbsint
#  
#  Created by Alexander Rudy on 2013-02-12.
#  Copyright 2012 Alexander Rudy. All rights reserved.
# 

try:
    import numpy.distutils
except ImportError:
    print("NUMPY>=1.6 is required for this package.")

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration(None,parent_package,top_path)
    config.set_options(ignore_setup_xxx_py=True,
                       assume_default_configuration=True,
                       delegate_options_to_subpackages=True,
                       quiet=True)
    config.add_subpackage('bsint','lib/bsint')
    # print "v",config.get_version()
    return config

if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(
    author = "Alexander Rudy",
    author_email = "arrudy@ucsc.edu",
    configuration = configuration,
    )
    