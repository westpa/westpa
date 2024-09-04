'''
YAML-based configuration files for WESTPA
'''

import os
import warnings
import yaml

try:
    from yaml import CLoader as YLoader
except ImportError:
    # fall back on Python implementation
    from yaml import Loader as YLoader

from . import extloader

NotProvided = object()


class ConfigValueWarning(UserWarning):
    pass


def warn_dubious_config_entry(entry, value, expected_type=None, category=ConfigValueWarning, stacklevel=1):
    if expected_type:
        warnings.warn(
            'dubious configuration entry {}: {} (expected type {})'.format(entry, value, expected_type), category, stacklevel + 1
        )
    else:
        warnings.warn('dubious configuration entry {}: {}'.format(entry, value), category, stacklevel + 1)


def check_bool(value, action='warn'):
    '''Check that the given ``value`` is boolean in type. If not, either
    raise a warning (if ``action=='warn'``) or an exception (``action=='raise'``).
    '''
    if action not in ('warn', 'raise'):
        raise ValueError('invalid action {!r}'.format(action))

    if not isinstance(value, bool):
        if action == 'warn':
            warnings.warn(
                'dubious boolean value {!r}, will be treated as {!r}'.format(value, bool(value)),
                category=ConfigValueWarning,
                stacklevel=2,
            )
        elif action == 'raise':
            raise ValueError('dubious boolean value {!r}, would be treated as {!r}'.format(value, bool(value)))
    else:
        return value


class ConfigItemMissing(KeyError):
    def __init__(self, key, message=None):
        self.key = key
        if message is None:
            message = 'configuration item missing: {!r}'.format(key)
        super().__init__(message)


class ConfigItemTypeError(TypeError):
    def __init__(self, key, expected_type, message=None):
        self.key = key
        self.expected_type = expected_type
        if message is None:
            message = 'configuration item {!r} must have type {!r}'.format(key, expected_type)
        super().__init__(message)


class ConfigValueError(ValueError):
    def __init__(self, key, value, message=None):
        self.key = key
        self.value = value
        if message is None:
            message = 'bad value {!r} for configuration item {!r}'.format(key, value)
        super().__init__(message)


class YAMLConfig:
    preload_config_files = ['/etc/westpa/westrc', os.path.expanduser('~/.westrc')]

    def __init__(self):
        self._data = {}

        for source in self.preload_config_files:
            self.update_from_file(source, required=False)

    def __repr__(self):
        return repr(self._data)

    def update_from_file(self, file, required=True):
        if isinstance(file, str):
            try:
                file = open(file, 'rt')
            except IOError:
                if required:
                    raise
                else:
                    return

        self._data.update(yaml.load(file, Loader=YLoader))
        file.close()

    def _normalize_key(self, key):
        if isinstance(key, str):
            key = (key,)
        else:
            try:
                key = tuple(key)
            except TypeError:
                key = (key,)
        return key

    def _resolve_object_chain(self, key, last=None):
        if last is None:
            last = len(key)
        objects = [self._data[key[0]]]
        for subkey in key[1:last]:
            objects.append(objects[-1][subkey])
        return objects

    def __getitem__(self, key):
        key = self._normalize_key(key)
        return self._resolve_object_chain(key)[-1]

    def __setitem__(self, key, value):
        key = self._normalize_key(key)

        try:
            objchain = self._resolve_object_chain(key, -1)
        except KeyError:
            # creation of a new (possibly nested) entry
            val = self._data
            for keypart in key[:-1]:
                try:
                    val = val[keypart]
                except KeyError:
                    val[keypart] = {}
                    val = val[keypart]
            try:
                val = val[key[-1]]
            except KeyError:
                val[key[-1]] = value
        else:
            objchain[-1][key[-1]] = value

    def __delitem__(self, key):
        key = self._normalize_key(key)
        objchain = self._resolve_object_chain(key, -1)
        del objchain[-1][key[-1]]

    def __contains__(self, key):
        try:
            self[key]
        except (KeyError, TypeError):
            return False
        else:
            return True

    def require(self, key, type_=None):
        '''Ensure that a configuration item with the given ``key`` is present. If
        the optional ``type_`` is given, additionally require that the item has that
        type.'''

        try:
            item = self[key]
        except KeyError:
            raise ConfigItemMissing(key)

        if type_ is not None:
            if not isinstance(item, type_):
                raise ConfigItemTypeError(item, type_)
        return item

    def require_type_if_present(self, key, type_):
        '''Ensure that the configuration item with the given ``key`` has the
        given type.'''

        try:
            item = self[key]
        except KeyError:
            return
        else:
            if not isinstance(item, type_):
                raise ConfigItemTypeError(item, type_)

    def coerce_type_if_present(self, key, type_):
        try:
            item = self[key]
        except KeyError:
            return
        else:
            if type_ is bool and not isinstance(item, bool):
                warn_dubious_config_entry(key, item, bool)
            self[key] = type_(item)

    def get(self, key, default=None):
        try:
            return self[key]
        except KeyError:
            return default

    def get_typed(self, key, type_, default=NotProvided):
        try:
            item = self[key]
        except KeyError as ke:
            if default is not NotProvided:
                item = default
            else:
                raise ke

        # Warn about possibly bad boolean
        if type_ is bool and not isinstance(item, bool):
            warn_dubious_config_entry(key, item, bool)

        return type_(item)

    def get_path(self, key, default=NotProvided, expandvars=True, expanduser=True, realpath=True, abspath=True):
        try:
            path = self[key]
        except KeyError as ke:
            if default is not NotProvided:
                path = default
            else:
                raise ke

        if expandvars:
            path = os.path.expandvars(path)
        if expanduser:
            path = os.path.expanduser(path)
        if realpath:
            path = os.path.realpath(path)
        if abspath:
            path = os.path.abspath(path)

        return path

    def get_pathlist(self, key, default=NotProvided, sep=os.pathsep, expandvars=True, expanduser=True, realpath=True, abspath=True):
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

        if expandvars:
            items = list(map(os.path.expandvars, items))
        if expanduser:
            items = list(map(os.path.expanduser, items))
        if realpath:
            items = list(map(os.path.realpath, items))
        if abspath:
            items = list(map(os.path.abspath, items))

        return items

    def get_python_object(self, key, default=NotProvided, path=None):
        try:
            qualname = self[key]
        except KeyError as ke:
            if default is not NotProvided:
                return default
            else:
                raise ke

        return extloader.get_object(qualname, path)

    def get_choice(self, key, choices, default=NotProvided, value_transform=None):
        try:
            value = self[key]
        except KeyError:
            if default is not NotProvided:
                value = default
            else:
                raise

        choices = set(choices)
        if value_transform:
            value = value_transform(value)
        if value not in choices:
            raise ConfigValueError(
                key,
                value,
                message='bad value {!r} for configuration item {!r} (valid choices: {!r})'.format(
                    value, key, tuple(sorted(choices))
                ),
            )
        return value
