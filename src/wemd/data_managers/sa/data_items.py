__metaclass__ = type

class DBDataItem:
    @classmethod
    def named_create(cls, name, value):
        newobj = cls(name)
        newobj.value = value
        return newobj
    
    def __repr__(self):
        return '<DBDataItem %r->%r>' % (self.name, self.value)
    
    def __init__(self, name, pvalue=None, cvalue=None, bvalue=None):
        self.name = name
        self.pvalue = pvalue
        self.cvalue = cvalue
        self.bvalue = bvalue
        
    def _get_value(self):
        if self.pvalue is not None:
            return self.pvalue
        elif self.bvalue is not None:
            return self.bvalue
        elif self.cvalue is not None:
            return self.cvalue
        else:
            return None
        
    def _set_value(self, obj):
        self.pvalue = self.cvalue = self.bvalue = None
        if obj is None: 
            return
        elif isinstance(obj, basestring):
            self.cvalue = obj
        elif isinstance(obj, buffer):
            self.bvalue = obj
        else:
            self.pvalue = obj
            
    def _del_value(self):
        self._set_value(None)
            
    value = property(_get_value, _set_value, _del_value)

