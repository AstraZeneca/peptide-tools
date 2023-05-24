from enum import Enum


class PKaType(str, Enum):
    ACIDIC = "acid"
    BASIC = "base"


class PKaMethod(str, Enum):
    ACD = "acd"
    PKA_MATCHER = "pkamatcher"
