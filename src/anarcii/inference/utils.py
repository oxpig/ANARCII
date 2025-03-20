import string

from torch.nn.utils.rnn import pad_sequence
from torch.utils.data import DataLoader

# All upper case letters, then all upper case letters doubled, then a space.
alphabet = (
    # All upper case letters.
    list(string.ascii_uppercase)
    # All upper case letters, doubled.
    + [2 * letter for letter in string.ascii_uppercase]
)

# Allowed and forbidden IMGT instertion
full_cdrs = list(range(27, 39)) + list(range(56, 66)) + list(range(105, 118))
cdr_instertion_starts = [32, 60, 111]
forbidden_cdr_insertions = [x for x in full_cdrs if x not in cdr_instertion_starts]
allowed_non_cdr_instertions = [x for x in range(1, 129) if x not in full_cdrs]


def collate_fn(batch):
    """
    Custom collate function for DataLoader.

    Ensures that each batch is padded dynamically to the longest sequence in
    the batch.
    """
    return pad_sequence(batch, batch_first=True, padding_value=0)


def dataloader(batch_size, tokenised_seqs):
    """
    Returns a DataLoader that batches sequences dynamically.

    Parameters:
    - batch_size (int): Number of sequences per batch.
    - tokenised_seqs (list of tensors): Tokenized sequences.

    Returns:
    - DataLoader: Batches of shape [batch_size, max_seq_len].
    """
    return DataLoader(tokenised_seqs, batch_size=batch_size, collate_fn=collate_fn)


def build_inward_list(length: int, start_num: int, end_num: int):
    """
    IMGT numbering is such that insertions are numbered differently depending
    on where they are found.

    Outside of a CDR: 44A, 44B, 44C
    Inside of a CDR: 111A, 111B, 112C, 112B, 112A

    The lang model simply outputs insertions with an X label. This fxn converts
    to a number with an insertion label and provides a list of tuples
    [(number,  letter), ...] based on the specified length,  start number, and
    end number.

    Parameters:
    - length (int): The # of X tokens in the X run.
    - start_num (int): The number that preceded the X run.
    - end_num (int): The after the X run (used conditionally if in a loop).

    Returns:
    - list of tuples: Each tuple consists of a number (either start_num or
    end_num) and a corresponding letter representing the instertion.

    Behavior (conditional on loop):
    - If start_num is within a predefined set of start numbers (cdrs), the function
    creates a structured sequence where the first half of the list uses
    start_num, and the second half transitions to end_num in a mirrored pattern.
    - If in the forbidden set then a value error is raised.
    - If start_num is not in the cdr start or forbidden set, the function simply cycles
    through the alphabet, pairing each letter with start_num.

    """
    result = []
    if int(start_num) in cdr_instertion_starts:
        # Calculate midpoint
        midpoint = length // 2  # Find the middle index by floor division
        # 5 becomes 2 ensures that we start at the higher number if uneven
        for i in range(midpoint):
            result.append((int(start_num), alphabet[i % len(alphabet)]))
        for i in range(midpoint, length):
            if length % 2 != 0:  # odd
                result.append(
                    (int(end_num), alphabet[(midpoint * 2 - i) % len(alphabet)])
                )
            elif length % 2 == 0:  # even
                result.append(
                    (int(end_num), alphabet[(midpoint * 2 - (i + 1)) % len(alphabet)])
                )
        return result

    elif int(start_num) in forbidden_cdr_insertions:
        raise ValueError("Forbidden cdr insertion predicted.")

    elif int(start_num) in allowed_non_cdr_instertions:
        for i in range(length):
            result.append((int(start_num), alphabet[i % len(alphabet)]))
        return result

    else:
        raise ValueError("Error in converting predicted insertions labels.")
