import rich.json
from _typeshed import Incomplete as Incomplete
from typing import Literal

console: Incomplete

class OB16Parameters:
    data: Incomplete
    def __init__(self, resolution: Literal['h0', 'h1']) -> None: ...
    def __rich__(self) -> rich.json.JSON: ...
