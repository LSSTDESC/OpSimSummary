class fieldSimlib(object):
    
    def __init__(self, fieldsimlibstring):
        self.simlibString = fieldsimlibstring 
        self._tup = self.split_fieldsimlib()
        self.headerstr = self._tup[0]
        self.datastring = self._tup[1]
        self.footer = self._tup[2]
    
    @property
    def fieldID(self):
    
    @property
    def meta(self):
        pass
    
    @property
    def data(self):
        pass
    
    def validate(self):
        pass
    
