# FIXME:  Refactor for new numbered sequence data structure.
def legacy_output(dt, verbose):
    if verbose:
        print(
            "Converting to legacy format. Three separate lists. \n",
            "A list of numberings, a list of all alignment details (contains, id, "
            "chain and score), and an empty list for hit tables. \n",
        )

    numbering, alignment_details, hit_tables = [], [], []
    for key, value in dt.items():
        if value["numbering"]:
            numbering.append(
                [(value["numbering"], value["query_start"], value["query_end"])]
            )
        else:
            numbering.append(None)

        # Changes for Ody needed here.
        new_dict = {}
        new_dict["chain_type"] = value["chain_type"]
        new_dict["scheme"] = value["scheme"]
        new_dict["query_name"] = key
        new_dict["query_start"] = value["query_start"]
        new_dict["query_end"] = value["query_end"]

        alignment_details.append([new_dict])

        hit_tables.append(None)

    return numbering, alignment_details, hit_tables
