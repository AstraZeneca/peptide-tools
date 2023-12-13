import os
import logging


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
