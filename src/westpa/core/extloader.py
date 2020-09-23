import importlib
import logging
import sys

log = logging.getLogger(__name__)


def load_module(module_name, path=None):
    """Load and return the given module, recursively loading containing packages as necessary."""
    if module_name in sys.modules:
        log.debug('module %r already loaded' % module_name)
        return sys.modules[module_name]

    if path is None:
        return importlib.import_module(module_name)

    spec_components = list(reversed(module_name.split('.')))
    qname_components = []
    mod_chain = []
    while spec_components:
        next_component = spec_components.pop(-1)
        qname_components.append(next_component)

        try:
            parent = mod_chain[-1]
            path = parent.__path__
        except IndexError:
            parent = None

        qname = '.'.join(qname_components)

        if qname in sys.modules:
            log.debug("Skipping preloaded module {}".format(qname))
            module = sys.modules[qname]
        else:
            log.debug('find_spec({!r}, {!r})'.format(qname, path))

            spec = importlib.machinery.PathFinder().find_spec(qname, path)

            if spec is None:
                raise ImportError(f'No module named {qname}')

            module = importlib.util.module_from_spec(spec)

            if spec.name not in sys.modules:
                sys.modules[spec.name] = module

            spec.loader.exec_module(module)

            # Make the module appear in the parent module's namespace
            if parent:
                setattr(parent, next_component, module)

            log.debug('module %r loaded' % qname)

        mod_chain.append(module)

    return module


def get_object(object_name, path=None):
    """Attempt to load the given object, using additional path information if given."""

    try:
        (modspec, symbol) = object_name.rsplit('.', 1)
    except ValueError:
        # no period found
        raise ValueError("object_name name must be in the form 'module.symbol'")

    log.debug('attempting to load %r from %r' % (symbol, modspec))
    module = load_module(modspec, path)

    # This will raise AttributeError (as expected) if the symbol is not in the module
    return getattr(module, symbol)
