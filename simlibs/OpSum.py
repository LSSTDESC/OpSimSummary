#!/usr/bin/env python


class simlibFieldInfo(object):
    '''
    Class to encapsulate the data in a single field and methods associated with
    the field
    '''

    def __init__(self):
        self._summary = 'Not Constructed\n'
    @property
    def meta(self):
        pass
    @property
    def data(self):
        pass
    def validate(self):
        pass

    def summary(self):
        return self._summary
