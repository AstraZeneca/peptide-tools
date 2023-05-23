import json

from pichemist.config import AA_SMARTS_SET_FILEPATH


class AASmarts(object):
    """Dynamic AA SMARTS object."""
    def __init__(self, data):
        for k, v in data.items():
            setattr(self, k, v)


def create_global_object_from_standardised_json(filepath):
    """
    Reads a SMARTS set from standardised JSON and converts it
    into an object with attributes compatible with the legacy logic
    of the software.

    """
    with open(filepath) as f:
        data = json.load(f)
    return AASmarts(data["smarts"])


AA_SMARTS_SET = create_global_object_from_standardised_json(
    AA_SMARTS_SET_FILEPATH).__dict__
