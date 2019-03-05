import sys

def version_testing_builtin():
    """A simple function ensuring compatibility between python2 and python3 on the __builtin__ module."""
    if sys.version[0] == '2':
        package = "__builtin__"
    elif sys.version[0] == '3':
        package = "builtins"
    else:
        raise Exception("This is not the python we're looking for (version {})".format(sys.version[0]))
    return package
