import os, sys, re, string
import logging
log = logging.getLogger('wemd.util.config_dict')

class ConfigError(EnvironmentError, ValueError):
    def __init__(self, message = '', exc = None):
        self.message = message
        self.exc = exc

    def __str__(self):
        if self.message:
            if self.exc:
                return '%s: %s' % (self.message, self.exc)
            else:
                return self.message
        elif self.exc:
            return str(self.exc)
        else:
            return ValueError.__str__(self)    

class ConfigDict(dict):
    true_values = set(('1', 'true', 't', 'yes', 'y'))
    false_values = set(('0', 'false', 'f', 'no', 'n'))
    default_list_split = re.compile(r'\s*,?\s*')

    defaults = {}
    
    def __init__(self, *args, **kwargs):
        self.update(self.defaults)
        super(ConfigDict,self).__init__(*args, **kwargs)
            
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
        for section in cparser.sections():
            log.debug('processing section %r' % section)
            for (key, value) in cparser.items(section):
                self['%s.%s' % (section, key)] = value
                log.debug('%s.%s = %s' % (section, key, value))
                
    def require(self, key):
        (section, name) = key.split('.', 1)
        if key not in self:
            raise ConfigError(("entry '%s' in section '%s' is required in "
                              +"configuration file %r")
                              % (name, section, self['__file__']))
    
    def require_all(self, keys):
        for key in keys:
            (section, name) = key.split('.', 1)
            if key not in self:
                raise ConfigError(("entry '%s' in section '%s' is required in "
                                   +"configuration file %r")
                                   % (name, section, self['__file__']))
                        
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
    
    
    def get_list(self, key, *args, **kwargs):
        if len(args) > 1:
            raise TypeError('unexpected positional argument encountered')
        
        type_ = kwargs.get('type')

        split = re.compile(kwargs.get('split', self.default_list_split))
        try:
            if type_:
                try:
                    return [type_(item) for item in split.split(self[key])]
                except ValueError, e:
                    raise ConfigError("%s in entry '%s'" % (e, key))
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
            
        if os.path.isabs(path):
            return os.path.normpath(path)
        else:
            return os.path.normpath(os.path.join(self['__dirname__'], path))
        
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
                    raise ConfigError('gzip compression not supported')
                
                return gzip.open(path, mode)
            elif ext == '.bz2':
                try:
                    import bz2
                except ImportError:
                    raise ConfigError('bzip2 compression not supported')
                
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
            raise ConfigError('invalid substitution template %r in %r: %s'
                              % (template.template, key, e))
        else:
            return template
        
    def get_python_object(self, key, *args):
        if len(args) > 1:
            raise TypeError('unexpected positional argument encountered')
        
        try:
            qualname = self[key]
        except KeyError, ke:
            try:
                return args[0]
            except IndexError:
                raise ke
        
        (module_name, symbol) = qualname.rsplit('.', 1)

        try:
            module = sys.modules[module_name]
        except KeyError:
            module = __import__(module_name)
            
        return getattr(module, symbol)