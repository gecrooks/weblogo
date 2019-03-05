
# TODO: replace with direct calls to pkg_resources
from weblogo.utils import resource_string, resource_stream, resource_filename


def data_string(name):
    return resource_string(__name__, 'data/' + name, __file__)


def data_stream(name):
    return resource_stream(__name__, 'data/' + name, __file__)


def data_filename(name):
    return resource_filename(__name__, 'data/' + name, __file__)
