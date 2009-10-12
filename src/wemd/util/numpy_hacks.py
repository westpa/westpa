class NumpyCmpSafeDict(dict):
    def __eq__(self, other):
        mykeys = set(self)
        try:
            otherkeys = set(other)
        except TypeError:
            return False
        
        if mykeys != otherkeys: return False
        
        for key in mykeys:
            myval = self[key]
            otherval = other[key]
            try:
                if myval == otherval:
                    continue
                else:
                    return False
            except ValueError, e:
                if 'ambiguous' in str(e):
                    try:
                        if (myval == otherval).all():
                            continue
                        else:
                            return False
                    except AttributeError:
                        raise e
