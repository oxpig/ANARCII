import pathlib

type Input = str | tuple[str, str] | list[str | tuple[str, str]] | dict[str, str]


def file_input(path: pathlib.Path) -> dict[str, str]: ...


def coerce_input(input_data: Input) -> dict[str, str]:
    """
    Coerce varied input sequence data formats into a dictionary.

    Accepts one or more peptide sequence strings, packaged up in a variety of ways,
    producing a dictionary with the sequences as values and the accompanying labels as
    keys:
    * A file path, which will be read to extract the sequence strings and their labels;
    * A single sequence string, which will be labelled with the key `sequence`;
    * A list of sequence strings which will be labelled sequentially with keys
      `sequence-1`, `sequence-2`, etc.;
    * A dictionary of name-sequence pairs, which will be returned unmodified;
    * A tuple of name-sequence pairs, or a list thereof, will be converted into a
      dictionary of the same.

    Args:
        input_data (Input): Peptide sequence strings, optionally labelled with names.

    Raises:
        TypeError: An unrecognised type of input data was provided.

    Returns:
        dict[str, str]: The values are sequence strings and the keys are names, which
                        are generated if not provided.
    """
    try:
        # Capture the cases list[tuple[str, str]] | dict[str, str],
        # containing name-sequence pairs.
        return dict(input_data)
    except ValueError:
        if isinstance(input_data, str):
            # Capture the case of file input.
            path = pathlib.Path(input_data)
            if path.suffix:
                return file_input(path)
            # Capture the case of a single peptide sequence string.
            return {"sequence": input_data}
        if isinstance(input_data, tuple):
            # Capture the case of a single name-sequence tuple[str, str].
            name, sequence = input_data
            return {name: sequence}
        if isinstance(input_data, list):
            # Capture the case of a list of peptide sequence strings, labelling
            # sequentially with `sequence-1`, `sequence-2`, etc..
            width = len(str(len(input_data)))
            return {
                f"sequence-{i:0{width}d}": seq for i, seq in enumerate(input_data, 1)
            }

    raise TypeError("Invalid input type.")
