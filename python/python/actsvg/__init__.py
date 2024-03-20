from .pyactsvg import *
from . import pyactsvg
from . import actsvg_json as json
from actsvg import io

setattr(io, "json", json)
