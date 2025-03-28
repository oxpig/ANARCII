from __future__ import annotations

import uuid
from collections.abc import Iterator
from itertools import chain, pairwise
from pathlib import Path
from typing import Any

import msgpack


def to_msgpack(object: Any, path: Path | str | None) -> Path:
    """
    Serialise an object to a file in MessagePack format.

    If no path is provided, a random UUID is generated and the filename is
    `anarcii-<UUID>.msgpack`.

    Args:
        object:  Any serialisable object.
        path:    Path to the output file.  If `None`, a filename is generated.

    Returns:
        Path to the output file.
    """
    path = path or f"anarcii-{uuid.uuid4()}.msgpack"
    with open(path, "wb") as f:
        msgpack.pack(object, f)

    return Path(path).absolute()


def from_msgpack_map(
    path: Path | str, chunk_size: int = 100 * 1024
) -> Iterator[dict[Any, Any]]:
    """
    Unpack a MessagePack map from a file.

    Yield a dictionary containing no more than `chunk_size` entries at a time.

    Args:
        path:        A file containing a MessagePack map as the first entry.
        chunk_size:  Maximum number of entries to yield at a time.

    Yields:
        A dictionary containing up to `chunk_size` entries at a time.
    """
    with open(path, "rb") as f:
        unpacker = msgpack.Unpacker(f, use_list=False)
        map_length = unpacker.read_map_header()
        for bounds in pairwise(chain(range(0, map_length, chunk_size), (map_length,))):
            yield {unpacker.unpack(): unpacker.unpack() for _ in range(*bounds)}
