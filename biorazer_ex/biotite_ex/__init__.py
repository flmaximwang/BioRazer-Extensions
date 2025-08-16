try:
    import biotite
except ImportError:
    raise ImportError(
        "This package requires biotite to be installed. "
        "Please install it using 'pip install biotite'."
    )
from . import structure
