import json
import re

import pandas as pd


def write_csv(ls, filename=None):
    # Works on the simple output
    rows = []
    for sublist in ls:
        row_dict = {}
        # Two rows for each sublist:
        # one for the first and
        # one for the second element of the tuples

        nums = sublist[0]
        name = sublist[1]["query_name"]
        chain = sublist[1]["chain_type"]
        score = sublist[1]["score"]

        row_dict["Name"] = name
        # Gives an F for the first letter in Fail.
        row_dict["Chain"] = chain
        row_dict["Score"] = score

        if chain == "F":
            rows.append(row_dict)
            continue

        for res in nums:
            key = str(res[0][0]) + res[0][1].strip()
            value = res[1].strip()
            row_dict[key] = value

        rows.append(row_dict)

    df = pd.DataFrame(rows)

    # Function to split into numeric and alphabetical parts
    def split_key(column):
        match = re.match(r"(\d+)([A-Z]?)", column)
        if match:
            num_part = int(match.group(1))  # Numeric part as integer
            letter_part = match.group(2)  # Alphabetical part as string
            return (num_part, letter_part)
        else:
            # For columns that don't match, return a tuple that will sort them
            # after the ones that do match.
            # Here, we use a large number and the column itself as a fallback.
            return (float("inf"), column)

    # # Sorting the columns list based on the custom key
    sorted_columns = sorted(
        [col for col in df.columns if col not in {"Name", "Chain", "Score"}],
        key=split_key,
    )
    sorted_columns = ["Name", "Chain", "Score"] + sorted_columns

    df = df[sorted_columns]

    if filename:
        df.to_csv(filename, index=False, header=True)

    return df


def write_json(data, file_path):
    with open(file_path, "w") as json_file:
        json.dump(data, json_file, indent=2)


def return_imgt_regions(ls):
    # Define IMGT regions with ranges
    imgt_regions = {
        "fw1": range(1, 27),
        "cdr1": range(27, 39),
        "fw2": range(39, 56),
        "cdr2": range(56, 66),
        "fw3": range(66, 105),
        "cdr3": range(105, 118),
        "fw4": range(118, 129),
    }

    ls_of_region_dicts = []

    for sublist in ls:
        nums = sublist[0]
        regions = {}

        # Populate regions by checking which IMGT region each number falls into
        for region_name, region_range in imgt_regions.items():
            regions[region_name] = [num for num in nums if num[0][0] in region_range]

        ls_of_region_dicts.append(regions)

    return ls_of_region_dicts
