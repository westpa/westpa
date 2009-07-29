class lazy_load(object):
    def __init__(self, backing_name):
        self.backing_name = backing_name
        
    def __get__(self, obj, objtype):
        try:
            return getattr(obj, self.backing_name)
        except AttributeError:
            val = obj.__lazy_load__(self.backing_name)
            setattr(obj, self.backing_name, val)
            return val
        
    def __set__(self, obj, val):
        setattr(obj, self.backing_name, val)
        
    def __delete__(self, obj):
        delattr(obj, self.backing_name)
