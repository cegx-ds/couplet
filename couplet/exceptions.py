#
# Copyright (C) 2022 Cambridge Epigenetix. All rights reserved.
#
# Define user-defined exceptions
class InsufficientPythonVersionError(Exception):
    """Raised when the Python version is insufficient."""

    def __init__(self, min_version, version):
        min_version_str = ".".join([str(i) for i in min_version])
        version_str = ".".join([str(i) for i in version])
        self.message = (
            "Python version (" + version_str + ") is insufficient. "
            "Required version: " + min_version_str + "."
        )
        super().__init__(self.message)


class IncompatibleErrorCodes(Exception):
    """Raised when trying to merge mismatch stats from data based on different rules."""

    def __init__(self, error_code_set):
        self.message = (
            "Attempting to merge data based on different rules. "
            "A total of "
            + str(len(error_code_set))
            + " sets of error codes were observed ("
            "where a single set was expected). The sets were: "
            + str(error_code_set)
            + "."
        )
        super().__init__(self.message)
