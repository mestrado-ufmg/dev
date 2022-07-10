from math import pi
from typing import Dict, Any, List, TypeVar, Callable, Type, cast

from pybird.models.enums import TailShape

deg2RadConst = pi / 180
rad2DegConst = 1 / deg2RadConst

T = TypeVar("T")

def from_dict(f: Callable[[Any], T], x: Any) -> Dict[str, T]:
    assert isinstance(x, dict)
    return { k: f(v) for (k, v) in x.items() }


def from_float(x: Any) -> float:
    assert isinstance(x, (float, int)) and not isinstance(x, bool)
    return float(x)


def to_float(x: Any) -> float:
    assert isinstance(x, float)
    return x


def from_str(x: Any) -> str:
    assert isinstance(x, str)
    return x

def from_list(f: Callable[[Any], T], x: Any) -> List[T]:
    assert isinstance(x, list)
    return [f(y) for y in x]


def to_class(c: Type[T], x: Any) -> dict:
    assert isinstance(x, c)
    return cast(Any, x).to_dict()

def to_tail_shape(x: Any) -> TailShape:
    assert isinstance(x, int)
    assert x <= 4

    if x == 1:
        return TailShape.rounded
    elif x == 2:
        return TailShape.square
    elif x == 3:
        return TailShape.pointed
    else:
        return TailShape.v

def from_tail_shape(c: Type[T], x: Any) -> int:

    assert isinstance(x, c)

    if x == TailShape.rounded:
        return 1
    elif x == TailShape.square:
        return 2
    elif x == TailShape.pointed:
        return 3
    else:
        return 4