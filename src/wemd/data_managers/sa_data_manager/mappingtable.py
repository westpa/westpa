import sqlalchemy
from sqlalchemy import insert, select, bindparam
from sqlalchemy.exceptions import OperationalError

class DictTableObject(object):
    def __init__(self, key, value):
        self.key = key
        self.value = value
    def __repr__(self):
        return repr( (self.key, self.value) )

class DictTableIface(object):
    def __init__(self, bind, table, key_field, value_field):
        self.bind = bind
        self.table = table
        self.key_field = key_field
        self.value_field = value_field
        
        self._get_item_query = select([self.table]).where(self.key_field == bindparam('keyval'))
        self._set_item_query = self.table.insert()
        self._del_item_query = self.table.delete(self.key_field == bindparam('keyval'))
        self._count_key_query = select([sqlalchemy.func.count(self.key_field)])
        self._get_items_query = select([self.table])
        self._get_keys_query = select([self.key_field])
        self._get_values_query = select([self.value_field])
        
    def __len__(self):
        return self.bind.execute(self._count_key_query).scalar()
    
    def __contains__(self, key):
        q = self._count_key_query.where(self.key_field == key)
        return bool(self.bind.execute(q).scalar() > 0)
    
    def __getitem__(self, key):
        rsl = self.bind.execute(self._get_item_query, {'keyval': key})
        row = rsl.fetchone()
        if row:
            return row[1]
        else:
            raise KeyError(key)
    
    def __setitem__(self, key, value):
        del self[key]
        self.bind.execute(self._set_item_query.values({self.key_field: key, 
                                                       self.value_field: value}))
    
    def __delitem__(self, key):
        self.bind.execute(self._del_item_query, {'keyval': key})
    
    def iterkeys(self):
        rsl = self.bind.execute(self._get_keys_query)
        row = rsl.fetchone()
        while row:
            yield row[0]
            row = rsl.fetchone()
        rsl.close()
        
    def itervalues(self):
        rsl = self.bind.execute(self._get_values_query)
        row = rsl.fetchone()
        while row:
            yield row[0]
            row = rsl.fetchone()
        rsl.close()
        
    def iteritems(self):
        rsl = self.bind.execute(self._get_items_query)
        row = rsl.fetchone()
        while row:
            yield (row[0], row[1])
            row = rsl.fetchone()
        rsl.close()
        
    def keys(self):
        return list(self.iterkeys())
    
    def values(self):
        return list(self.itervalues())
    
    def items(self):
        return list(self.iteritems())

    def get(self, key, default=None):
        try:
            return self[key]
        except KeyError:
            return default
    
    def setdefault(self, key, value):
        try:
            return self[key]
        except KeyError:
            self[key] = value
            return value
        
    def update(self, *args, **kwargs):
        updatedict = dict(*args)
        updatedict.update(**kwargs)
        for (k,v) in updatedict.iteritems():
            self[k] = v
