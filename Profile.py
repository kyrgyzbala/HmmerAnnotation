'''
Created on Jun 3, 2014

@author: sanjarbek
'''

import params as p
import os


class HmmProfile:
    def __init__(self, hmm_file_path, profile_type=None, order=[], name='not-named'):
        self.file_path = hmm_file_path
        self.source = os.path.splitext(os.path.basename(hmm_file_path))[0]
        self.raw_lines = open(hmm_file_path).readlines()
        self.profile_type=profile_type
        self.order=order
        self.name=name

