__metaclass__ = type
from sqlalchemy.orm.collections import MappedCollection

class DBDataItem:
    def __repr__(self):
        return '<%s %r->%r>' % (self.__class__.__name__, self.name, self.value)
    
    def __init__(self, name, value = None, 
                 pvalue=None, cvalue=None, bvalue=None):
        self.name = name
        if value is not None:
            self.value = value
        else:
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
    
class KeyFunc(object):
    def __init__(self, attr):
        self.attr = attr
    def __call__(self, obj):
        return getattr(obj, self.attr)
    
class _dict(dict): pass
        
class DBDataDict(_dict, MappedCollection):
    def __init__(self, *args, **kwargs):
        MappedCollection.__init__(self, keyfunc = KeyFunc('name'))
        dict.__init__(self, *args, **kwargs)
