from __future__ import annotations

import uuid
from collections.abc import Iterator
from itertools import chain, pairwise
from pathlib import Path
from typing import Any, BinaryIO

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


def _open_msgpack_map_file(
    f: BinaryIO, chunk_size: int = 100 * 1024
) -> Iterator[dict[Any, Any]]:
    """
    Unpack a MessagePack map from a file object opened in binary mode.

    Args:
        f:           A file object containing a MessagePack map as the first entry.
        chunk_size:  Maximum number of entries to yield at a time.

    Yields:
        A dictionary containing no more than `chunk_size` entries at a time.
    """
    unpacker = msgpack.Unpacker(f, use_list=False)
    map_length = unpacker.read_map_header()
    for bounds in pairwise(chain(range(0, map_length, chunk_size), (map_length,))):
        yield {unpacker.unpack(): unpacker.unpack() for _ in range(*bounds)}


def from_msgpack_map(
    path: Path | str, chunk_size: int = 100 * 1024
) -> Iterator[dict[Any, Any]]:
    """
    Unpack a MessagePack map from a file.

    Args:
        path:        A file containing a MessagePack map as the first entry.
        chunk_size:  Maximum number of entries to yield at a time.

    Yields:
        A dictionary containing no more than `chunk_size` entries at a time.
    """
    with open(path, "rb") as f:
        yield from _open_msgpack_map_file(f, chunk_size)
