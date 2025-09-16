import json
import os
import sys
from io import StringIO


script_dir = os.path.dirname(os.path.realpath(__file__))
examples_dir = os.path.join(script_dir, "examples")


class TestError(Exception):
    """Raised when need to break the tests."""

    __test__ = False


def stdout_to_variable(func, *args, **kwargs):
    """
    Assign sys.stdout somewhere else and plug the
    result variable to that pipe to collect the print.

    """
    tmp_stdout = sys.stdout
    sys.stdout = result = StringIO()
    func(*args, **kwargs)
    sys.stdout = tmp_stdout
    return result.getvalue()


def jsonify_output(output):
    return json.loads(json.dumps(output))
