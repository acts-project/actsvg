from .pyactsvg import *
from . import pyactsvg
from . import actsvg_json as json
from . import actsvg_uproot as uproot
from actsvg import io

setattr(io, "json", json)

