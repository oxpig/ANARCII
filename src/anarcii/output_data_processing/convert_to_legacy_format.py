# FIXME:  Refactor for new numbered sequence data structure.
def convert_output(ls, format, verbose):
    if format == "simple":
        return ls

    elif format == "legacy":
        if verbose:
            print(
                "Converting to legacy format. Three separate lists. \n",
                "A list of numberings, a list of all alignment details (contains, id, "
                "chain and score), and an empty list for hit tables. \n",
            )

        numbering, alignment_details, hit_tables = [], [], []
        for x in ls:
            if x[0]:
                numbering.append([(x[0], x[1]["query_start"], x[1]["query_end"])])
            else:
                numbering.append(None)

            # Changes for Ody needed here.
            new_dict = x[1]
            new_dict["species"] = None
            new_dict["scheme"] = "imgt"

            alignment_details.append([new_dict])

            hit_tables.append(None)

        return numbering, alignment_details, hit_tables
