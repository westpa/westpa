import os, sys, re, string
import logging
import extloader
from collections import OrderedDict
from datetime import timedelta
log = logging.getLogger('wemd.util.config_dict')

NotProvided = object()

class ConfigItemMissing(KeyError):
    def __init__(self, key, message=None):
        self.missing_key = key
        if message is None:
            message = 'configuration item missing: {!r}'.format(key)
        super(ConfigItemMissing,self).__init__(message)

class ConfigDict(OrderedDict):
    true_values = set(('1', 'true', 't', 'yes', 'y'))
    false_values = set(('0', 'false', 'f', 'no', 'n'))
    default_list_split = re.compile(r'\s*,?\s*')
    
    re_interval_float_unit = re.compile(r'((?:\d+\.)?\d+)\s+(second|minute|hour|day|sec|min|hr|s|m|h|d)s?')
    re_interval_dhms = re.compile(r'(?:\d+:)?\d+:\d\d:\d\d')

    defaults = {}
    
    def __init__(self, *args, **kwargs):
        self.update(self.defaults)
        super(ConfigDict,self).__init__(*args, **kwargs)
        
    def update_from_object(self, obj, prefix='args'):
        udict = {'{}.{}'.format(prefix,key): val for (key,val) in obj.__dict__.viewitems() if not key.startswith('_')}
        if log.isEnabledFor(logging.DEBUG):
            log.debug('updating with {!r}'.format(udict))
        self.update(udict)
            
    def read_config_file(self, config_filename):
        from ConfigParser import SafeConfigParser
                
        log.debug('opening %r' % config_filename)
        config_file = open(config_filename, 'r')
        self.filename = self['__file__'] = os.path.abspath(config_filename)
        self.dirname  = self['__dir__'] = os.path.dirname(self.filename)

        cparser = SafeConfigParser({'__file__': self.filename,
                                    '__dir__': self.dirname})
        cparser.optionxform = str
        cparser.readfp(config_file)
        
        
        defaults = cparser.defaults()
        log.debug('DEFAULT section')
        for (key, value) in defaults.viewitems():
            log.debug('%s = %r' % (key, value))
        default_keys = set(defaults.viewkeys())
        
        for section in cparser.sections():
            log.debug('reading section %r' % section)
            for (key, value) in cparser.items(section):
                self['%s.%s' % (section, key)] = value
                if key in default_keys and defaults[key] == value:
                    log.debug('%s.%s propagated from [DEFAULT]' % (section,key))
                else:
                    log.debug('%s.%s = %r' % (section, key, value))
                
    def require(self, key):
        try:
            return self[key]
        except KeyError:
            raise ConfigItemMissing(key)
    
    def require_all(self, keys):
        return map(self.require, keys)
                        
    def _get_typed(self, key, type_, *args):
        
        try:
            return type_(self[key])
        except KeyError, ke:
            try:
                return args[0]
            except IndexError:
                return ke

    def get_int(self, key, *args):
        if len(args) > 1:
            raise TypeError('unexpected positional argument encountered')

        try:
            return int(self[key])
        except KeyError, ke:
            try:
                return args[0]
            except IndexError:
                raise ke
    
    def get_bool(self, key, *args, **kwargs):
        if len(args) > 1:
            raise TypeError('unexpected positional argument encountered')

        try:
            strict = kwargs['strict']
            del kwargs['strict']
        except KeyError:
            strict = True
            
        if kwargs:
            raise TypeError('unexpected keyword argument %r encountered'
                            % kwargs.pop()[0] )
            
        try:
            val = str(self[key]).lower()
        except KeyError, ke:
            try:
                return args[0]
            except IndexError:
                raise ke

        if val in self.false_values:
            return False
        else:
            if strict:
                if val in self.true_values:
                    return True
                else:
                    raise ValueError('invalid boolean literal %r' % val)
            else:
                return True

    def get_float(self, key, *args):
        if len(args) > 1:
            raise TypeError('unexpected positional argument encountered')

        try:
            return float(self[key])
        except KeyError, ke:
            try:
                return args[0]
            except IndexError:
                raise ke
    
    def get_interval(self, key, default=NotProvided, type=None):        
        try:
            inttxt = self[key].lower()
        except KeyError as ke:
            if default is not NotProvided:
                inttxt = default
            else:
                raise ke
        
        if not inttxt:
            return None
        
        match = self.re_interval_float_unit.match(inttxt) or self.re_interval_dhms.match(inttxt)
        if not match:
            ValueError('invalid time interval %r' % inttxt)
        groups = match.groups()
        
        multipliers = {'s': 1, 'm': 60, 'h': 3600, 'd': 86400}
        if groups:
            # matched re_interval_float_unit, returning (quantity, unit)
            # since unit is uniquely determined by its first character, just examine that and
            # multiply appropriately
            seconds = float(groups[0]) * multipliers[groups[1][0]]
        else:
            qtys = map(float, match.group(0).split(':'))
            seconds = qtys[-1] + 60*qtys[-2] + 3600*qtys[-3]
            try:
                seconds += qtys[-4] * 86400
            except IndexError:
                pass
            
        if type is None:
            return timedelta(seconds=seconds)
        else:
            return type(seconds)        
    
    def get_list(self, key, *args, **kwargs):
        if len(args) > 1:
            raise TypeError('unexpected positional argument encountered')
        
        type_ = kwargs.get('type')

        split = re.compile(kwargs.get('split', self.default_list_split))
        try:
            if type_:
                [type_(item) for item in split.split(self[key])]
            else:
                return split.split(self[key])
        except KeyError, ke:
            try:
                return args[0]
            except IndexError:
                raise ke
            
    def get_path(self, key, *args):
        if len(args) > 1:
            raise TypeError('unexpected positional argument encountered')

        try:
            path = self[key]
        except KeyError, ke:
            try:
                path = args[0]
            except IndexError:
                raise ke
        
        path = os.path.expanduser(os.path.expandvars(path))    
        if os.path.isabs(path):
            return os.path.normpath(path)
        else:
            return os.path.normpath(os.path.join(self['__dir__'], path))
        
    def get_pathlist(self, key, default=NotProvided, sep=os.pathsep, 
                     expandvars = True, expanduser = True, realpath = True, abspath = True):
        try:
            paths = self[key]
        except KeyError as ke:
            if default is not NotProvided:
                paths = default
            else:
                raise ke
        
        try:
            items = paths.split(sep)
        except AttributeError:
            # Default must have been something we can't process, like a list or None
            # Just pass it through, since enforcing a restriction on what kind of
            # default is passed is probably more counterproductive than any poor programming
            # practice it encourages.
            return paths
        
        if expandvars: items = map(os.path.expandvars, items)
        if expanduser: items = map(os.path.expanduser, items)
        if realpath:   items = map(os.path.realpath, items)
        if abspath:    items = map(os.path.abspath, items)
        
        return items
                
    def get_file_object(self, key, default_path=None, mode='rb', 
                        compression_allowed=True):
        try:
            path = self[key]
        except KeyError, ke:
            if default_path is None:
                raise ke
            else:
                path = default_path
        
        if compression_allowed:
            (bn, ext) = os.path.splitext(path)
            if ext == '.gz':
                try:
                    import gzip
                except ImportError:
                    raise ValueError('gzip compression not supported')
                
                return gzip.open(path, mode)
            elif ext == '.bz2':
                try:
                    import bz2
                except ImportError:
                    raise ValueError('bzip2 compression not supported')
                
                return bz2.BZ2File(path, mode[0]+'U')
        return open(path, mode)
        
    def get_compiled_template(self, key, *args):
        if len(args) > 1:
            raise TypeError('unexpected positional argument encountered')
        
        try:
            template = self[key]
        except KeyError, ke:
            try:
                template = args[0]
            except IndexError:
                raise ke
                
        if not template: 
            return None
        elif not isinstance(template, string.Template):
            template = string.Template(template)
        
        try:
            template.safe_substitute(dict())
        except ValueError, e:
            raise ValueError('invalid substitution template %r in %r: %s'
                              % (template.template, key, e))
        else:
            return template
        
    def get_python_object(self, key, default=NotProvided, path=None):
        try:
            qualname = self[key]
        except KeyError as ke:
            if default is not NotProvided:
                return default
            else:
                raise ke
        
        return extloader.get_object(qualname, path)
    
    def get_python_callable(self, key, default=NotProvided, path=None):
        try:
            qualname = self[key]
        except KeyError as ke:
            if default is not NotProvided:
                return default
            else:
                raise ke
            
        fn = extloader.get_object(qualname, path)
        if not callable(fn):
            raise TypeError('%s is not callable' % qualname)

        
