import argparse
import logging
import os


def get_logger(name):
    """
    Creates a logger where level is
    determined by the env LOGGING_LEVEL.

    """
    log = logging.getLogger(name)
    if os.environ.get("LOGGING_LEVEL") == "DEBUG":
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    return log


def str2bool(v):
    """
    Converts a string to a boolean.
    """
    if isinstance(v, bool):
        return v
    if v.lower() in ("yes", "true", "t", "y", "1"):
        return True
    elif v.lower() in ("no", "false", "f", "n", "0"):
        return False
    else:
        raise argparse.ArgumentTypeError("Boolean value expected.")
