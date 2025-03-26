from .schemes_utils import conversion_function


def convert_number_scheme(numbered_seqs_dict, scheme):
    """Renumber a dict of IMGT seqs with new scheme.

    This takes a dict of IMGT numbered sequences.
    It works out if each sequence is a heavy or light chain
    Defines the scheme to be applied
    Then calls the conversion function on that sequence
    """

    converted_seqs = {}
    for nm, dt in numbered_seqs_dict.items():
        chain_call = dt["chain_type"]
        chain = "heavy" if chain_call == "H" else "light"

        if scheme.lower() == "imgt":
            converted_seqs[nm] = conversion_function(dt, scheme.lower())
        else:
            scheme_name = scheme.lower() + "_" + chain
            converted_seqs[nm] = conversion_function(dt, scheme_name)

    return converted_seqs
