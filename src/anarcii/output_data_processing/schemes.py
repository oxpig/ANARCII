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
        if dt["numbering"] and dt["chain_type"] != "F":
            chain_call = dt["chain_type"]
            chain = "heavy" if chain_call == "H" else "light"

            if scheme.lower() == "imgt":
                converted_seqs[nm] = conversion_function(dt, scheme.lower())
            else:
                scheme_name = scheme.lower() + "_" + chain
                converted_seqs[nm] = conversion_function(dt, scheme_name)

        else:
            # Sequence is a fail - could not be numbered
            dt["scheme"] = scheme
            converted_seqs[nm] = dt

    return converted_seqs
